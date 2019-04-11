import re
import numpy as np

sequon_re = re.compile("(?=(N[^PX][ST]))")
sequon_re_modded = re.compile("(?=(N(?P<mod>\[[A-Za-z0-9\.+\-]+\])?[^PX][ST]))")


class Sequence:
    def __init__(self, seq, metadata=None):
        self.seq = seq
        if metadata:
            self.metadata = metadata
        else:
            self.metadata = dict()

    def find_with_regex(self, motif, ignore=None):
        pattern = motif
        new_str = ""
        if ignore is not None:
            for i in range(len(ignore)):
                if not ignore[i]:
                    new_str += self.seq[i]
        else:
            new_str = self.seq
        for i in pattern.finditer(new_str):
            for m in range(1, len(i.groups()) + 1):
                yield slice(i.start(m), i.end(m))

    def gaps(self):
        s = np.arange(len(self.seq))
        for i in range(len(s)):
            if self.seq[i] == '-':
                s[i] = True
            else:
                s[i] = False
        return s

    def count(self, char, start, end):
        return self.seq.count(char, start, end)


class Peptide:
    def __init__(self, seq):
        self.seq = seq
        self.modded = False
        self._prep()
        if not self.modded:
            self.stripped_sequon = self.seq
        self.glycocat = 0

    def _prep(self):
        if "[" in self.seq:
            self.modded = True
            self.stripped_sequon = re.sub("\[[A-Za-z0-9\.+\-]+\]", "", self.seq).replace('-', '')

    def map_seq(self, full_source, extra_len=2):
        seq_length = len(full_source)
        # stripped_sequon = re.sub("\[[A-Za-z0-9\.+\-]+\]", "", self.seq)
        original_position = full_source.find(self.stripped_sequon)

        stop_position = original_position + len(self.stripped_sequon)
        temp_seq = original_position + len(self.stripped_sequon) + extra_len
        if temp_seq <= seq_length:
            extra_seq = full_source[stop_position: len(self.stripped_sequon) + extra_len]
        else:
            extra_seq = full_source[stop_position: seq_length]
        return original_position, stop_position, extra_seq
