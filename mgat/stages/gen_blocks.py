import os
from pathos import multiprocessing
from mgat.utils.message import Message
from mgat.utils.fasta_io import FastaIO


def __split_fasta(in_fa, win_size, step_size, out_fa):
    Message.info("Splitting %s." % in_fa)
    fasta_io = FastaIO(in_fa)
    fasta_io.load()
    fasta_io.sliding_window_to_file(win_size, step_size, out_fa)


def gen_blocks(qry_dir, split_dir, win_size, step_size, threads):
    if not os.path.exists(split_dir):
        os.makedirs(split_dir)
    split_tasks = []
    split_files = []
    for qry_fn in os.listdir(qry_dir):
        in_fa = os.path.join(qry_dir, qry_fn)
        qry_name_list = qry_fn.split(".")[:-1]
        qry_name_suf = qry_fn.split(".")[-1]
        qry_name_list.append("split")
        qry_name_list.append(qry_name_suf)
        out_fa = os.path.join(split_dir, ".".join(qry_name_list))
        if os.path.exists(out_fa):
            Message.info("Split file %s already exists, skipping..." % out_fa)
            split_files.append(in_fa)
        else:
            split_tasks.append([in_fa, out_fa])

    if split_tasks:
        pool = multiprocessing.Pool(
            processes=threads if threads < len(split_tasks) else len(split_tasks)
        )
        res = []
        for in_fa, out_fa in split_tasks:
            r = pool.apply_async(
                __split_fasta,
                (
                    in_fa,
                    win_size,
                    step_size,
                    out_fa,
                ),
            )
            res.append([in_fa, r])
        pool.close()
        pool.join()

        for in_fa, r in res:
            try:
                r.get()
                split_files.append(in_fa)
            except Exception as e:
                Message.warn("Exception caught with {}: {}".format(in_fa, e))

    if len(split_files) == len(os.listdir(qry_dir)):
        return True
    else:
        return False
