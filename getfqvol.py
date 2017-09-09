import os
import sys
import gzip

class FqVol(object):
    def __init__(self, fq):
        self.fq = fq

    @staticmethod
    def __get_read_length(fq):
        with gzip.open(fq, "rb") as inf:
            fline = inf.readline()
            sline = inf.readline()
            orilen = len(sline)
            return orilen

    @staticmethod
    def __get_ori_num(fq):
        orirnumcmd = "zcat {fq} | wc -l | cut -d ' ' -f 1".format(fq=fq)
        orirnum = float(os.popen(orirnumcmd).read().strip())
        return orirnum

    def value(self):
        fq = self.fq

        orilen = FqVol.__get_read_length(fq)
        orinum = FqVol.__get_ori_num(fq)

        orivol = orilen * 2.0 * orinum / 4.0 / (10.0 ** 9.0)

        return orivol
