import sys
import os
import logging

class ComVcf(object):
    def __init__(self, vcf1, vcf2):
        self.vcf1 = vcf1
        self.vcf2 = vcf2

    @staticmethod
    def __load_info(vcf):
        resdict = {}
        with open(vcf, "r") as invcf:
            for line in invcf:
                if not line.strip().startswith("#"):
                    infofields = line.strip().split()
                    chrid, pos, ids, ref, alt = [infofields[x] for x in [0, 1, 2, 3, 4]]
                    keys = "||".join([chrid, pos, ids, ref, alt])
                    resdict[keys] = line.strip()

    def comvcf(self):
        vcf1 = self.vcf1
        vcf2 = self.vcf2

        dict1 = ComVcf.__load_info(vcf1)
        dict2 = ComVcf.__load_info(vcf2)

        return cmp(dict1, dict2)

