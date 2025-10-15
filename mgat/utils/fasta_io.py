class FastaIO:
    def __init__(self, in_fa):
        self.__in_fa = in_fa
        self.fa_db = {}
        self.len_db = {}

    def load(self):
        with open(self.__in_fa) as fin:
            for line in fin:
                if line[0] == ">":
                    sid = line.strip().split()[0][1:]
                    self.fa_db[sid] = []
                else:
                    self.fa_db[sid].append(line.strip())

        for sid in self.fa_db:
            self.fa_db[sid] = "".join(self.fa_db[sid])
            self.len_db[sid] = len(self.fa_db[sid])

    def sliding_window_to_file(self, window_size, step_size, out_fa):
        with open(out_fa, "w") as fout:
            for sid in self.fa_db:
                for _ in range(0, len(self.fa_db[sid]) - window_size + 1, step_size):
                    sub_seq = self.fa_db[sid][_ : _ + window_size]
                    sub_sid = "%s::%d::%d" % (
                        sid,
                        _ + 1,
                        (
                            _ + window_size
                            if _ + window_size < self.len_db[sid]
                            else self.len_db[sid]
                        ),
                    )
                    fout.write(">%s\n%s\n" % (sub_sid, sub_seq))
