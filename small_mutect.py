import os
import sys
import argparse
import logging
from mutect import rushver, oriver
from setlog import setlog
from comvcf import ComVcf

def opt():
    args = argparse.ArgumentParser(description="test small data for mutect")
    args.add_argument("-a", "--anadir", help="The absolute path of analysis directory", dest="anadir",
                      required=True)
    return args.parse_args()

class Mutect(object):
    def __init__(self, anadir):
        self.bamdir = os.path.join(anadir, "bams")
        self.beddir = os.path.join(anadir, "rois")
        self.outdir = os.path.join(anadir, "output")

        self.run()

    @staticmethod
    def __get_bam(bamdir):
        casedict = {}
        controldict = {}

        for p, ds, fs in os.walk(bamdir):
            for f in fs:
                if f.endswith(".bam"):
                    if "cancer" in f:
                        sample = f.rstrip("_cancer.bam")
                        casedict[sample] = os.path.join(p, f)
                        controldict[sample] = os.path.join(p, "{}_normal.bam".format(sample))
        return casedict, controldict

    @staticmethod
    def __get_bed(beddir, sample):
        bedpath = os.path.join(beddir, "{}.bed".format(sample))
        if not os.path.isfile(bedpath):
            logging.error("{} is not existed, please check!".format(bedpath))
            bedpath = ""
        return bedpath

    def run(self):
        bamdir = self.bamdir
        beddir = self.beddir
        outdir = self.outdir
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        casedict, controldict = Mutect.__get_bam(bamdir)

        for sample in casedict.keys():
            if sample in controldict.keys():
                logging.info("Treating sample: {}".format(sample))
                subdir = os.path.join(outdir, sample)
                if not os.path.isdir(subdir):
                    os.makedirs(subdir)

                casebam = casedict[sample]
                controlbam = controldict[sample]
                bed = Mutect.__get_bed(beddir, sample)
                if not os.path.isfile(os.path.join(subdir, "time.log")):
                    rushtimestr, rushvcf = rushver(case=casebam, control=controlbam, anadir=subdir, targetbed=bed)
                    oritimestr, orivcf = oriver(case=casebam, control=controlbam, anadir=subdir, targetbed=bed)
                    with open(os.path.join(subdir, "time.log"), "w") as ouf:
                        ouf.write("\n".join([rushtimestr, oritimestr]))

                    comvcfres = ComVcf(vcf1=rushvcf, vcf2=orivcf).comvcf()
                    if comvcfres == 0:
                        logging.info("Both of vcf files is common for {}".format(sample))
                    else:
                        logging.error("Both of vcf files is diff for {}".format(sample))
                else:
                    logging.debug("{} existed, skipping".format(os.path.join(subdir, "time.log")))

def main():
    setlog()
    args = opt()
    anadir = args.anadir

    Mutect(anadir=anadir)

if __name__ == '__main__':
    main()

