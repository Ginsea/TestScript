import os
import sys
import argparse
import logging
import time
import multiprocessing

from mutect import rushver, oriver, getbam
from setlog import setlog
from comvcf import ComVcf

def opt():
    args = argparse.ArgumentParser(description="run multiple beds")
    args.add_argument("-a", "--anadir", dest="anadir", required=True, help="The current analysis directory")
    return args.parse_args()

def main():
    args = opt()
    anadir = args.anadir

    setlog()

    bamdir = os.path.join(anadir, "bam")
    beddir = os.path.join(anadir, "bed")

    casebam, controlbam = getbam(bamdir)

    with open(os.path.join(anadir, "{}.log".format(time.strftime("%H%M%S", time.localtime(time.time())))), "a") as logf:
        for p, ds, fs in os.walk(beddir):
            for f in fs:
                if f.endswith(".bed"):
                    pool = multiprocessing.Pool(processes=2)
                    reslist = []
                    bedfile = os.path.join(p, f)
                    reslist.append(pool.apply_async(rushver, args=(casebam, controlbam, anadir, bedfile)))
                    reslist.append(pool.apply_async(oriver, args=(casebam, controlbam, anadir, bedfile)))
                    pool.close()
                    pool.join()

                    vcflist = []
                    timestrlist = []
                    for subres in reslist:
                        timestr, vcf = subres.get()
                        vcflist.append(vcf)
                        timestrlist.append(timestr)

                    vcf1, vcf2 = vcflist
                    comvcf = ComVcf(vcf1, vcf2).comvcf()
                    logf.write("{}\t{}\t{}\n".format(os.path.basename(bedfile), "\t".join(timestrlist), comvcf))

                    
if __name__ == '__main__':
    main()
