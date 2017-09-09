import os
import sys
import argparse
import multiprocessing
import random
import logging
import subprocess

from setlog import setlog

def findbam(anadir):
    logging.info("searching bam in {}".format(anadir))
    reslist = []
    for p, ds, fs in os.walk(anadir):
        for f in fs:
            if f.endswith("_cancer_sort_markdup_realign_recal_ori.bam"):
                reslist.append(os.path.join(p, f))
    return reslist

class TransData(object):
    def __init__(self, anadir):
        self.anadir = anadir

    @staticmethod
    def __get_bai(bam):
        bamdir = os.path.dirname(bam)
        bamname = os.path.splitext(os.path.basename(bam))[0]
        baifile = os.path.join(bamdir, "{}.bai".format(bamname))
        if not os.path.isfile(baifile):
            logging.error("{} is not existed, please check!".format(baifile))
            exit(1)
        return baifile

    @staticmethod
    def __get_norbam(casebam):
        casedir = os.path.dirname(casebam)
        normaldir = os.path.join(casedir, "../../normal/5_recal_bam/")
        samplename = os.path.basename(casebam).rstrip("_cancer_sort_markdup_realign_recal_ori.bam")
        normalbam = os.path.join(normaldir, "{}_normal_sort_markdup_realign_recal_ori.bam".format(samplename))
        if not os.path.isfile(normalbam):
            logging.error("{} is not existed, please check!")
            exit(1)
        return normalbam


    def transcmd(self):
        anadir = self.anadir

        casebamlist = findbam(anadir)
        randomcasebam = random.sample(casebamlist, 50)

        transcmdlist = []

        for casebam in randomcasebam:
            normalbam = TransData.__get_norbam(casebam)
            casebai = TransData.__get_bai(casebam)
            normalbai = TransData.__get_bai(normalbam)

            for subfile in [casebam, normalbam, casebai, normalbai]:
                cmd = "scp {} chenzhx@10.0.23.100:/data/chenzhx/FPGA/anadir/2017-09-07/clinic/bams && " \
                      "touch {}.SUCCESS".format(subfile, os.path.basename(subfile))
                transcmdlist.append(cmd)

        return transcmdlist

def runcmd(cmd):
    subp = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    logging.info("Transferring {}".format(cmd.split()[1]))
    logging.debug("[cmd] {}".format(cmd))
    subp.wait()
    if subp.returncode != 0:
        logging.error("Error reported in {} transferring".format(cmd.split()[1]))
    else:
        logging.info("Transferring {} finished".format(cmd.split()[1]))

def opt():
    args = argparse.ArgumentParser(description="transfer clinical bam file from 10.0.1.101 to 10.0.23.100")
    args.add_argument("-a", "--anadir", help="The absolute path of analysis directory", dest="anadir", required=True)
    return args.parse_args()

def main():
    args = opt()

    anadir = args.anadir
    setlog()

    cmdlist = TransData(anadir).transcmd()

    pool = multiprocessing.Pool(processes=25)

    for subcmd in cmdlist:
        pool.apply_async(runcmd, args=(subcmd, ))

    pool.close()
    pool.join()

if __name__ == '__main__':
    main()