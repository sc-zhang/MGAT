import os
from mgat.utils.runner import Minimap2Runner
from mgat.utils.message import Message


def mapping_with_minimap2(ref_dir, qry_dir, mapping_dir, minimap_params, threads):
    if not os.path.exists(mapping_dir):
        os.makedirs(mapping_dir)
    mapping_tasks = []
    mapped_files = []
    for ref_fn in os.listdir(ref_dir):
        full_ref_fn = os.path.join(ref_dir, ref_fn)
        for qry_fn in os.listdir(qry_dir):
            full_qry_fn = os.path.join(qry_dir, qry_fn)
            bam_fn = "%s_vs_%s.bam" % (
                ".".join(ref_fn.split(".")[:-1]),
                ".".join(qry_fn.split(".")[:-1]),
            )
            full_bam_fn = os.path.join(mapping_dir, bam_fn)
            mapped_files.append(full_bam_fn)
            if os.path.exists(full_bam_fn):
                Message.info("Bam file %s exists, skipping..." % full_bam_fn)
            else:
                mapping_tasks.append([full_ref_fn, full_qry_fn, full_bam_fn])

    if mapping_tasks:
        mapper = Minimap2Runner()
        for ref_fn, qry_fn, bam_fn in mapping_tasks:
            mapper.mapping(ref_fn, qry_fn, bam_fn, threads, minimap_params)

    for bam_fn in mapped_files:
        if not os.path.exists(bam_fn):
            return False
    return True
