import os
from pathos import multiprocessing
from mgat.utils.message import Message


def __classify(sub_iden_db, threshold):
    classify_db = {}
    for chrn in sub_iden_db:
        classify_db[chrn] = {}
        for sp in sub_iden_db[chrn]:
            classify_db[chrn][sp] = []
            iden_list = []
            if sub_iden_db[chrn][sp]:
                for ref_name in sub_iden_db[chrn][sp]:
                    iden_list.append([ref_name, sub_iden_db[chrn][sp][ref_name]])

                if len(iden_list) == 1:
                    best_ref_name = iden_list[0][0]
                    best_iden = iden_list[0][1]
                    classify_db[chrn][sp] = [best_ref_name, best_iden]
                else:
                    iden_list = sorted(iden_list, key=lambda x: x[1], reverse=True)
                    best_ref_name = iden_list[0][0]
                    best_iden = iden_list[0][1]

                    if best_iden < iden_list[1][1] * threshold:
                        classify_db[chrn][sp] = ["Undetermined", best_iden]
                    else:
                        classify_db[chrn][sp] = [best_ref_name, best_iden]
            else:
                classify_db[chrn][sp] = ["Undetermined", 0]
    return classify_db


def classify_blocks(iden_db, classify_dir, threshold, threads):
    if not os.path.exists(classify_dir):
        os.makedirs(classify_dir)

    classify_tasks = []
    for qry_name in iden_db:
        classify_fn = os.path.join(classify_dir, qry_name + ".cla")
        if os.path.exists(classify_fn):
            Message.info("Classify of %s already exists, skipping" % qry_name)
        else:
            classify_tasks.append(qry_name)

    classify_db = {}
    if classify_tasks:
        pool = multiprocessing.Pool(
            processes=threads if threads < len(classify_tasks) else len(classify_tasks)
        )
        res = []
        for qry_name in classify_tasks:
            r = pool.apply_async(
                __classify,
                (
                    iden_db[qry_name],
                    threshold,
                ),
            )
            res.append([qry_name, r])
        pool.close()
        pool.join()

        for qry_name, r in res:
            try:
                classify_db[qry_name] = r.get()
            except Exception as e:
                Message.warn(
                    "Exception caught when classify {}: {}".format(qry_name, e)
                )

        if classify_db:
            for qry_name in classify_db:
                classify_fn = os.path.join(classify_dir, qry_name + ".cla")
                Message.info("Writing classify file %s" % classify_fn)
                with open(classify_fn, "w") as fout:
                    fout.write("#Chrom\tStart\tClass\tIdentity\n")
                    for chrn in sorted(classify_db[qry_name]):
                        for sp in sorted(classify_db[qry_name][chrn]):
                            ref_name, iden = classify_db[qry_name][chrn][sp]
                            fout.write("%s\t%d\t%s\t%f\n" % (chrn, sp, ref_name, iden))

    for fn in os.listdir(classify_dir):
        qry_name = ".".join(fn.split(".")[:-1])
        if qry_name in classify_db:
            continue
        classify_db[qry_name] = {}
        classify_fn = os.path.join(classify_dir, fn)
        Message.info("Loading exists classify file %s" % classify_fn)
        with open(classify_fn, "r") as fin:
            for line in fin:
                if line[0] == "#":
                    continue
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[1])
                cla = data[2]
                iden = float(data[3])
                if chrn not in classify_db[qry_name]:
                    classify_db[qry_name][chrn] = {}
                classify_db[qry_name][chrn][sp] = [cla, iden]
    return classify_db
