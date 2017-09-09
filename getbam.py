# -*- coding:utf-8 -*-

import os
import sys
import re
import subprocess
import multiprocessing
import logging
import argparse

from varibales import *
from setlog import setlog

def run_cmd(cmd):
    logging.info("[cmd] {}".format(cmd))
    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subp.wait()
    if subp.returncode != 0:
        logging.error("{} failed, please check!".format(cmd))
        exit(1)
    else:
        logging.debug("{} finished".format(cmd))

class Magicdict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class RunBwa(object):
    def __init__(self, anadir):
        self.anadir = anadir

    @staticmethod
    def __get_lane(instr):
        regstr = re.compile(r"[\d]+_[A-Za-z0-9]+_[A-Za-z0-9]+_(L[\d]+)_[A-Za-z0-9]+\-[\d]+\-v[A-Za-z0-9]+")
        lane = regstr.search(instr).groups()[0]
        return lane

    @staticmethod
    def __get_fq(fqdir):
        fqdict = Magicdict()
        for p, ds, fs in os.walk(fqdir):
            for f in fs:
                if f.endswith("_1.fq.gz"):
                    lane = RunBwa.__get_lane(str(f))
                    if "case" in str(p) or "cancer" in str(p):
                        fqdict["case"][lane]["fq1"] = os.path.join(p, f)
                    elif "normal" in str(p) or "control" in str(p):
                        fqdict["normal"][lane]["fq1"] = os.path.join(p, f)
                elif f.endswith("_2.fq.gz"):
                    lane = RunBwa.__get_lane(str(f))
                    if "case" in str(p) or "cancer" in str(p):
                        fqdict["case"][lane]["fq2"] = os.path.join(p, f)
                    elif "normal" in str(p) or "control" in str(p):
                        fqdict["normal"][lane]["fq2"] = os.path.join(p, f)
        return fqdict

    @staticmethod
    def __run_bwa(fqdict, anadir):
        cmdtemp = "{bwa} mem -t 6 -R '@RG\\tID:{label}\\tSM:{label}\\tLB:{label}Lib\\tPU:runname\\tCN:GenePlus\\tPL:illumina' " \
                  "{ref} {fq1} {fq2} > {outdir}/{outsam} 2> {logdir}/{log} && " \
                  "{samtools} view -@ 6 -bS {outdir}/{outsam} > {outdir}/{outbam} && " \
                  "{samtools} sort -@ 6 -m 1536M {outdir}/{outbam} -o {outdir}/{sortbam} && " \
                  "touch {logdir}/{slog}"

        outdir = os.path.join(anadir, "output/align")
        logdir = os.path.join(anadir, "log/align")
        for subdir in [outdir, logdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        pool = multiprocessing.Pool(processes=6)
        for keys in fqdict.keys():
            for lane in fqdict[keys].keys():
                sam = "{}_{}.sam".format(keys, lane)
                bam = "{}_{}.bam".format(keys, lane)
                sortbam = "{}_{}_sort.bam".format(keys, lane)
                log = "{}_{}.bwa.log".format(keys, lane)
                slog = "{}_{}.bwa.SUCCESS".format(keys, lane)
                fq1 = fqdict[keys][lane]["fq1"]
                fq2 = fqdict[keys][lane]["fq2"]
                cmd = cmdtemp.format(bwa=bwa, label=keys, ref=genome, fq1=fq1, fq2=fq2,
                                     outdir=outdir, outsam=sam, logdir=logdir, log=log,
                                     slog=slog, samtools=samtools, outbam=bam, sortbam=sortbam)
                pool.apply_async(run_cmd, args=(cmd, ))
        pool.close()
        pool.join()

        return outdir

    @staticmethod
    def __merge_bam(bamdir):
        casebam = []
        controlbam = []
        for p, ds, fs in os.walk(bamdir):
            for f in fs:
                if f.endswith("_sort.bam"):
                    if f.startswith("case") or f.startswith("cancer"):
                        casebam.append(os.path.join(p, f))
                    elif f.startswith("normal") or f.startswith("control"):
                        controlbam.append(os.path.join(p, f))

        casemergebam = os.path.join(bamdir, "case_sort.bam")
        controlmergebam = os.path.join(bamdir, "control_sort.bam")

        cmd = "{samtools} merge -@ 6 -f {mergebam} {allbam} && " \
              "{samtools} index {mergebam}"

        casemergecmd = cmd.format(samtools=samtools, mergebam=casemergebam,
                                  allbam=" ".join(casebam))
        controlmergecmd = cmd.format(samtools=samtools, mergebam=controlmergebam,
                                     allbam=" ".join(controlbam))

        logging.info("[cmd] {}".format(casemergecmd))
        logging.info("[cmd] {}".format(controlmergecmd))

        casesubp = subprocess.Popen(casemergecmd, shell=True, stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE)
        controlsubp = subprocess.Popen(controlmergebam, shell=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        casesubp.wait()
        controlsubp.wait()

        if casesubp.returncode != 0:
            logging.error("Error for cmd: {}, info is {}".format(casemergecmd, casesubp.stderr.read()))
        else:
            logging.debug("cmd: {} has been finished!.".format(casemergecmd))

        if controlsubp.returncode != 0:
            logging.error("Error for cmd: {}, info is {}".format(controlmergecmd, controlsubp.stderr.read()))
        else:
            logging.debug("cmd: {} has been finished!.".format(controlmergecmd))

        return casemergebam, controlmergebam

    def getbam(self):
        anadir = self.anadir

        fqdir = os.path.join(anadir, "fqpath")

        fqdict = RunBwa.__get_fq(fqdir)
        bamdir = RunBwa.__run_bwa(fqdict=fqdict, anadir=anadir)
        casebam, controlbam = RunBwa.__merge_bam(bamdir)

        return casebam, controlbam


def runpicard(bam, anadir):
    outdir = os.path.join(anadir, "output/picard")
    logdir = os.path.join(anadir, "log/picard")
    tempdir = os.path.join(anadir, "temp/picard")
    for subdir in [outdir, logdir, tempdir]:
        if not os.path.isdir(subdir):
            os.makedirs(subdir)

    outbam = os.path.join(outdir, "{}_markdup.bam".format(os.path.splitext(os.path.basename(bam))[0]))
    outmetrics = os.path.join(outdir, "{}_markdup.metrics".format(os.path.splitext(os.path.basename(bam))[0]))
    flog = os.path.join(logdir, "picard.log")
    slog = os.path.join(logdir, "picard.SUCCESS")

    cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms8g  -Xmx16g  " \
          "-Djava.io.tmpdir={temp} -jar {picard} MarkDuplicates I={inbam} O={outbam} " \
          "METRICS_FILE={outmetrics} ASO=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 " \
          "VALIDATION_STRINGENCY=LENIENT > {logfile} 2>&1 && " \
          "touch {slog} && " \
          "{samtools} index {outbam}".format(java=java, temp=tempdir, picard=picard,
                                             inbam=bam, outbam=outbam, outmetrics=outmetrics,
                                             logfile=flog, slog=slog, samtools=samtools,
                                             )

    run_cmd(cmd)

def opt():
    args = argparse.ArgumentParser(description="get bam file")
    args.add_argument("-a", "--anadir", help="The absolute path of current analysis work folder",
                      required=True, dest="anadir")
    args.add_argument("-b", "--bam", help="bam file path, default is None", default=None, dest="bam")
    return args.parse_args()

def main():
    args = opt()
    setlog()
    anadir = args.anadir
    bam = args.bam

    if not bam:
        casebam, controlbam = RunBwa(anadir=anadir).getbam()

        for subbam in [casebam, controlbam]:
            runpicard(bam=subbam, anadir=anadir)
    else:
        runpicard(bam=bam, anadir=anadir)

if __name__ == '__main__':
    main()









        






