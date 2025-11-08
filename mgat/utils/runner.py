import os.path
from subprocess import Popen, PIPE
from mgat.utils.message import Message


class Runner:

    def __init__(self):
        self._cmd = ""
        self._res = ""
        self._err = ""

    def set_command(self, cmd):
        self._cmd = '"%s"' % cmd
        self._res = ""
        self._err = ""

    def print_command(self):
        Message.info(self._cmd)

    def run(self):
        p = Popen(self._cmd, stdout=PIPE, stderr=PIPE, shell=True, encoding="utf-8")
        self._res, self._err = p.communicate()

    def get_result(self):
        return self._res

    def get_err(self):
        return self._err


class Minimap2Runner(Runner):
    def mapping(
        self, ref_file, qry_file, out_bam, threads, minimap_params="-ax asm5 --eqx"
    ):
        if os.path.getsize(ref_file) >= 4e9:
            if "--split-prefix" not in minimap_params:
                minimap_params += " --split-prefix %s" % (out_bam + ".minimap.idx")
        self._cmd = (
            "minimap2 {} -t {} {} {} | samtools sort -@ {} -o {} - -T {}.tmp".format(
                minimap_params,
                threads,
                ref_file,
                qry_file,
                threads,
                out_bam,
                out_bam,
            )
        )
        self.print_command()
        self.run()
        with open(out_bam + ".mapping.log", "w") as fout:
            fout.write("%s\n%s\n" % (self.get_result(), self.get_err()))
        Message.info("Minimap2 finished")

        Message.info("Indexing")
        self._cmd = "samtools index -@ %d %s" % (threads, out_bam)
        self.print_command()
        self.run()
        with open(out_bam + ".index.log", "w") as fout:
            fout.write("%s\n%s\n" % (self.get_result(), self.get_err()))
        Message.info("Indexing finished")
