import os
import sys
import argparse
import logging

from mutect import rushver
from setlog import setlog

def opt():
    args = argparse.ArgumentParser(description="test rush mutect through huge data")
    args.add_argument("-b", "--bams", help="bam directory", dest="bams", required=True)
    args.add_argument("-r", "--region", help="target region", dest="region", required=True)
    return args.parse_args()

def getbam(bamdir):
    normaldict = {}
    casedict = {}
    for p, ds, fs in os.walk(bamdir):
        for f in fs:
            if f.endswith("_normal_sort_markdup_realign_recal_ori.bam"):
                samplename = f.rstrip("_normal_sort_markdup_realign_recal_ori.bam")
                normaldict[samplename] = os.path.join(p, f)
            elif f.endswith("_cancer_sort_markdup_realign_recal_ori.bam"):
                samplename = f.rstrip("_cancer_sort_markdup_realign_recal_ori.bam")
                casedict[samplename] = os.path.join(p, f)
    return casedict, normaldict

def main():
    args = opt()
    bamdir = args.bamdir
    region = args.region

    logfile = "command.log"
    setlog(logfile=logfile)

    anadir = os.getcwd()

    resfile = os.path.join(anadir, "hugetest.res")

    casedict, normaldict = getbam(bamdir)

    for samplename in casedict.keys():
        if samplename in normaldict.keys():
            casebam = casedict[samplename]
            normalbam = normaldict[samplename]
            logging.info("case bam:{}\tnormal bam:{}".format(casebam, normalbam))
            timestr, vcf = rushver(case=casebam, control=normalbam, anadir=anadir, targetbed=region)
            with open(resfile, "a") as rf:
                rf.write("{}\t{}\t{}\t{}\n".format(samplename, os.path.getsize(casebam),
                                                   os.path.getsize(normalbam),
                                                   timestr))

if __name__ == '__main__':
    main()

