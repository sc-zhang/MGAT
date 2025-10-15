import argparse
import os
import shutil
import sys
from mgat.utils.message import Message
from mgat.stages.gen_blocks import gen_blocks
from mgat.stages.mapping import mapping_with_minimap2
from mgat.stages.identity import identify_blocks_ancestor
from mgat.stages.classify import classify_blocks
from mgat.stages.chain_classify import chain_classify_blocks


def get_opts():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r", "--reference", help="Directory of reference genomes", required=True
    )
    parser.add_argument(
        "-q", "--query", help="Directory of query genomes", required=True
    )
    parser.add_argument(
        "-w", "--window", help="Window size, default=5000", default=5000, type=float
    )
    parser.add_argument(
        "-s", "--step", help="Step size, default=5000", default=5000, type=float
    )
    parser.add_argument(
        "-o", "--output", help="Directory of output files", required=True
    )
    parser.add_argument(
        "--minimap_args",
        help='Arguments for minimap, default="-ax asm10 --eqx"',
        default="-ax asm10 --eqx",
    )
    parser.add_argument(
        "--threshold",
        help="Threshold for classifying blocks by identity, "
        "means best identity should higher than second best identity * threshold, "
        "default=1.5",
        type=float,
        default=1.5,
    )
    parser.add_argument(
        "--chain",
        help="Block counts for classifying undetermined blocks, the blocks with this counts before current block and"
        "the blocks with this counts after current block would be used, default=5",
        type=int,
        default=5,
    )
    parser.add_argument(
        "-t", "--threads", help="Number of threads, default=10", default=10, type=int
    )
    return parser.parse_args()


def pipeline(opts):
    ref_dir = os.path.abspath(opts.reference)
    qry_dir = os.path.abspath(opts.query)
    out_dir = os.path.abspath(opts.output)
    win_size = opts.window
    step_size = opts.step
    minimap_args = opts.minimap_args
    threshold = opts.threshold
    chain = opts.chain
    threads = opts.threads

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    Message.info("Step1. Splitting query genomes...")
    split_dir = os.path.join(out_dir, "01.split")
    status = gen_blocks(qry_dir, split_dir, win_size, step_size, threads)
    if status:
        Message.info("Splitting finished")
    else:
        Message.error("Splitting failed")
        sys.exit(1)

    Message.info("Step2. Running minimap2...")
    mapping_dir = os.path.join(out_dir, "02.mapping")
    status = mapping_with_minimap2(
        ref_dir, split_dir, mapping_dir, minimap_args, threads
    )
    if status:
        Message.info("Mapping finished")
    else:
        Message.error("Mapping failed")
        sys.exit(1)

    Message.info("Step3. Getting blocks identities...")
    iden_dir = os.path.join(out_dir, "03.identities")
    block_iden_db = identify_blocks_ancestor(mapping_dir, iden_dir, threads)
    if block_iden_db:
        Message.info("Getting blocks identities finished")
    else:
        Message.error("Getting blocks identities failed")
        sys.exit(-1)

    Message.info("Step4. Classifying blocks...")
    classify_dir = os.path.join(out_dir, "04.classify")
    classify_db = classify_blocks(block_iden_db, classify_dir, threshold, threads)
    if classify_db:
        Message.info("Classify blocks finished")
    else:
        Message.error("Classify blocks failed")
        sys.exit(-1)

    Message.info("Step5. Classifying undetermined blocks...")
    chain_classify_dir = os.path.join(out_dir, "05.chain_classify")
    status = chain_classify_blocks(
        classify_db, chain_classify_dir, win_size, step_size, threshold, chain, threads
    )
    if status:
        Message.info("Chain classify blocks finished")
    else:
        Message.error("Chain classify blocks failed")
        sys.exit(-1)

    Message.info("Step6. Generating final results...")
    result_dir = os.path.join(out_dir, "06.results")
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    for fn in os.listdir(chain_classify_dir):
        src = os.path.join(chain_classify_dir, fn)
        dst = os.path.join(result_dir, fn)
        shutil.copy(src, dst)
    Message.info("Generating final results finished")

    Message.info("Finished")


def main():
    opts = get_opts()
    pipeline(opts)
