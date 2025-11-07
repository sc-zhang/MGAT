import os
import pysam
from pathos import multiprocessing
from mgat.utils.message import Message


def __calc_identity(bam_file, ref_chrn):
    Message.info("Getting identity with %s on %s" % (bam_file, ref_chrn))
    iden_db = {}
    with pysam.AlignmentFile(bam_file, "rb") as fin:
        for line in fin.fetch(contig=ref_chrn):
            qry_id = line.query_name
            chrn, sp, ep = qry_id.split("::")
            sp = int(sp)
            ep = int(ep)
            block_size = ep - sp + 1
            match_size = 0
            for aln_type, base_cnt in line.cigartuples:
                if aln_type == 0 or aln_type == 7:
                    match_size += base_cnt

            iden = match_size * 1.0 / block_size
            if chrn not in iden_db:
                iden_db[chrn] = {}
            if sp not in iden_db[chrn] or iden > iden_db[chrn][sp]:
                iden_db[chrn][sp] = iden
    return iden_db


def identify_blocks_ancestor(split_dir, mapping_dir, iden_dir, threads):
    if not os.path.exists(iden_dir):
        os.makedirs(iden_dir)

    iden_tasks = []
    for fn in os.listdir(mapping_dir):
        if fn.endswith(".bam"):
            tmp = ".".join(fn.split(".")[:-1])
            ref_name, qry_name = tmp.split("_vs_")
            bam_fn = os.path.join(mapping_dir, fn)
            iden_fn = os.path.join(iden_dir, qry_name + ".iden")
            if os.path.exists(iden_fn):
                Message.info(
                    "Identity of %s vs %s already exists, skipping"
                    % (qry_name, ref_name)
                )
                continue
            with pysam.AlignmentFile(bam_fn, "rb") as fin:
                ref_chrs = fin.references
            for chrn in ref_chrs:
                iden_tasks.append([qry_name, ref_name, bam_fn, chrn])

    merged_db = {}
    if iden_tasks:
        pool = multiprocessing.Pool(
            processes=threads if threads < len(iden_tasks) else len(iden_tasks)
        )
        res = []
        for qry_name, ref_name, bam_fn, chrn in iden_tasks:
            r = pool.apply_async(
                __calc_identity,
                (
                    bam_fn,
                    chrn,
                ),
            )
            res.append([qry_name, ref_name, bam_fn, chrn, r])
        pool.close()
        pool.join()

        for qry_name, ref_name, bam_fn, ref_chrn, r in res:
            try:
                if qry_name not in merged_db:
                    merged_db[qry_name] = {}
                iden_db = r.get()
                for chrn in iden_db:
                    if chrn not in merged_db[qry_name]:
                        merged_db[qry_name][chrn] = {}
                    for sp in iden_db[chrn]:
                        if sp not in merged_db[qry_name][chrn]:
                            merged_db[qry_name][chrn][sp] = {}
                        merged_db[qry_name][chrn][sp][ref_name] = iden_db[chrn][sp]

            except Exception as e:
                Message.warn(
                    "Exception caught with {} on {}: {}".format(bam_fn, ref_chrn, e)
                )

        if merged_db:
            for qry_name in merged_db:
                iden_fn = os.path.join(iden_dir, qry_name + ".iden")
                Message.info("Writing identify file %s" % iden_fn)
                with open(iden_fn, "w") as fout:
                    fout.write("#Chrom\tStart\tRef\tIdentity\n")
                    for chrn in sorted(merged_db[qry_name]):
                        for sp in sorted(merged_db[qry_name][chrn]):
                            for ref_name in sorted(merged_db[qry_name][chrn][sp]):
                                fout.write(
                                    "%s\t%d\t%s\t%f\n"
                                    % (
                                        chrn,
                                        sp,
                                        ref_name,
                                        merged_db[qry_name][chrn][sp][ref_name],
                                    )
                                )

    for fn in os.listdir(iden_dir):
        qry_name = ".".join(fn.split(".")[:-1])
        if qry_name in merged_db:
            continue
        merged_db[qry_name] = {}
        iden_fn = os.path.join(iden_dir, fn)
        Message.info("Loading exist identity %s" % iden_fn)
        with open(iden_fn, "r") as fin:
            for line in fin:
                if line[0] == "#":
                    continue
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[1])
                ref_name = data[2]
                iden = float(data[3])
                if chrn not in merged_db[qry_name]:
                    merged_db[qry_name][chrn] = {}
                if sp not in merged_db[qry_name][chrn]:
                    merged_db[qry_name][chrn][sp] = {}
                if ref_name != "None":
                    merged_db[qry_name][chrn][sp][ref_name] = iden

    # Fill merged_db with blocks which could not be mapped
    for fn in os.listdir(split_dir):
        full_fn = os.path.join(split_dir, fn)
        qry_name = ".".join(fn.split(".")[:-1])
        if qry_name not in merged_db:
            merged_db[qry_name] = {}
        with open(full_fn, "r") as fin:
            for line in fin:
                if line.strip() and line[0] == ">":
                    chrn, sp, _ = line.strip()[1:].split("::")
                    sp = int(sp)
                    if chrn not in merged_db[qry_name]:
                        merged_db[qry_name][chrn] = {}
                    if sp not in merged_db[qry_name][chrn]:
                        merged_db[qry_name][chrn][sp] = {}
    return merged_db
