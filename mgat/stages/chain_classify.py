import os
from pathos import multiprocessing
from mgat.utils.message import Message


def __chain_classify(sub_classify_db, chain_block_counts):
    chain_classify_db = {}
    for chrn in sub_classify_db:
        chain_classify_db[chrn] = {}
        block_list = sorted(sub_classify_db[chrn])
        block_cnt = len(block_list)

        iter_cnt = 1
        change_cnt = 0
        while iter_cnt == 1 or change_cnt > 0:
            Message.info("%s, iteration %d" % (chrn, iter_cnt))
            iter_cnt += 1
            change_cnt = 0
            for i in range(block_cnt):
                if sub_classify_db[chrn][block_list[i]][0] == "Undetermined":
                    cnt_db = {}
                    for _ in range(
                        max(0, i - chain_block_counts),
                        min(block_cnt, i + chain_block_counts),
                    ):
                        if _ == i:
                            continue
                        block_class = sub_classify_db[chrn][block_list[_]][0]
                        if block_class == "Undetermined":
                            continue
                        if block_class not in cnt_db:
                            cnt_db[block_class] = 0
                        cnt_db[block_class] += 1

                    cnt_list = []
                    for block_class in cnt_db:
                        cnt_list.append([block_class, cnt_db[block_class]])
                    if cnt_list:
                        cnt_list = sorted(cnt_list, key=lambda x: x[1], reverse=True)
                        if len(cnt_list) == 1 or cnt_list[0][1] > cnt_list[1][1]:
                            sub_classify_db[chrn][block_list[i]][0] = cnt_list[0][0]
                            change_cnt += 1
            Message.info("%s classified %d undetermined blocks" % (chrn, change_cnt))

        for sp in sub_classify_db[chrn]:
            chain_classify_db[chrn][sp] = sub_classify_db[chrn][sp][0]
    return chain_classify_db


def chain_classify_blocks(
    classify_db,
    chain_classify_dir,
    win_size,
    step_size,
    threshold,
    chain_block_counts,
    threads,
):
    if not os.path.exists(chain_classify_dir):
        os.makedirs(chain_classify_dir)

    chain_classify_tasks = []
    chain_classify_result_set = set()
    for qry_name in classify_db:
        chain_classify_fn = os.path.join(chain_classify_dir, qry_name + ".chain.cla")
        chain_classify_result_set.add(chain_classify_fn)
        if os.path.exists(chain_classify_fn):
            Message.info("Chain classify of %s already exists, skipping" % qry_name)
        else:
            chain_classify_tasks.append(qry_name)

    chain_classify_db = {}
    if chain_classify_tasks:
        pool = multiprocessing.Pool(
            processes=(
                threads
                if threads < len(chain_classify_tasks)
                else len(chain_classify_tasks)
            )
        )
        res = []
        for qry_name in chain_classify_tasks:
            r = pool.apply_async(
                __chain_classify,
                (
                    classify_db[qry_name],
                    chain_block_counts,
                ),
            )
            res.append([qry_name, r])
        pool.close()
        pool.join()

        for qry_name, r in res:
            try:
                chain_classify_db[qry_name] = r.get()
            except Exception as e:
                Message.warn(
                    "Exception caught when chain classify {}: {}".format(qry_name, e)
                )

    if chain_classify_db:
        for qry_name in chain_classify_db:
            chain_classify_fn = os.path.join(
                chain_classify_dir, qry_name + ".chain.cla"
            )
            Message.info("Writing chain classify file %s" % chain_classify_fn)
            header = [
                "## Window size: %d" % win_size,
                "## Step size: %d" % step_size,
                "## Identity threshold: %f" % threshold,
                "## Chain block counts: %d" % chain_block_counts,
            ]
            with open(chain_classify_fn, "w") as fout:
                fout.write("%s\n" % ("\n".join(header)))
                fout.write("#Chrom\tStart\tClass\n")
                for chrn in sorted(chain_classify_db[qry_name]):
                    for sp in sorted(chain_classify_db[qry_name][chrn]):
                        fout.write(
                            "%s\t%d\t%s\n"
                            % (chrn, sp, chain_classify_db[qry_name][chrn][sp])
                        )

    for fn in chain_classify_result_set:
        if not os.path.exists(fn):
            return False
    return True
