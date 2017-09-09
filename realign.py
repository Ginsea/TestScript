import os
import sys
import subprocess
import time
import logging
import argparse

from varibales import *
from setlog import setlog

def runcmd(cmd, log=None):
    utime = 0
    stime = time.time()
    logging.info("[cmd] {}".format(cmd))
    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.debug("[pid] {}".format(subp.pid))
    subp.wait()
    if subp.returncode != 0:
        logging.error("Error for cmd: {}, please check {}".format(cmd, log))
        exit(1)
    else:
        etime = time.time()
        utime = (etime - stime) / 60
        logging.debug("cmd: {} has been finished, used time: {} minutes".format(cmd, utime))

    if log:
        with open(log, "a") as lf:
            lf.write(cmd)

    return utime

class RealignerTargetCreator(object):
    def __init__(self, anadir, label):
        self.outdir = anadir
        self.label = label

    @staticmethod
    def __get_bam(bamdir):
        cancer, normal = "", ""
        for p, ds, fs in os.walk(bamdir):
            for f in fs:
                if f.endswith("_markdup.bam"):
                    if "case" in f or "cancer" in f:
                        cancer = os.path.join(p, f)
                    elif "normal" in f or "control" in f:
                        normal = os.path.join(p, f)
        if not cancer or not normal:
            logging.error("can't find bam in {}".format(bamdir))
            exit(1)
        return cancer, normal

    def rushver(self):
        anadir = self.outdir
        bamdir = os.path.join(anadir, "output/picard")

        outdir = os.path.join(anadir, "output/realign/rush")
        logdir = os.path.join(anadir, "log/realign/rush")
        tempdir = os.path.join(anadir, "temp/realign/rush")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        casebam, normalbam = RealignerTargetCreator.__get_bam(bamdir)
        rtcout = os.path.join(outdir, "coRealign.intervals")
        rtclog = os.path.join(logdir, "RealignerTargetCreator.log")
        rtcslog = os.path.join(logdir, "RealignerTargetCreator.SUCCESS")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g  " \
              "-Djava.io.tmpdir={temp} -Djava.library.path={lib} -jar {gatk} -T RealignerTargetCreator " \
              "-I {normalbam} -I {casebam} -o {rtcout} -known {mills} -known {dbsnp} -L {targetbed} -R {genome} " \
              "-nt 12 -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment -dt NONE > {rtclog} 2>&1 && " \
              "touch {rtcslog}".format(java=java, temp=tempdir, lib=gatklib, gatk=newgatk, normalbam=normalbam,
                                       casebam=casebam, rtcout=rtcout, mills=mills, dbsnp=dbsnp, targetbed=targetbed,
                                       genome=genome, rtclog=rtclog, rtcslog=rtcslog)

        utime = runcmd(cmd=cmd, log=rtclog)

        return utime

    def oriver(self):
        anadir = self.outdir
        bamdir = os.path.join(anadir, "output/picard")

        outdir = os.path.join(anadir, "output/realign/ori")
        logdir = os.path.join(anadir, "log/realign/ori")
        tempdir = os.path.join(anadir, "temp/realign/ori")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        casebam, normalbam = RealignerTargetCreator.__get_bam(bamdir)
        rtcout = os.path.join(outdir, "coRealign.intervals")
        rtclog = os.path.join(logdir, "RealignerTargetCreator.log")
        rtcslog = os.path.join(logdir, "RealignerTargetCreator.SUCCESS")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g  " \
              "-Djava.io.tmpdir={temp} -jar {gatk} -T RealignerTargetCreator " \
              "-I {normalbam} -I {casebam} -o {rtcout} -known {mills} -known {dbsnp} -L {targetbed} -R {genome} " \
              "-nt 12 -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment -dt NONE > {rtclog} 2>&1 && " \
              "touch {rtcslog}".format(java=java, temp=tempdir, gatk=oldgatk, normalbam=normalbam,
                                       casebam=casebam, rtcout=rtcout, mills=mills, dbsnp=dbsnp, targetbed=targetbed,
                                       genome=genome, rtclog=rtclog, rtcslog=rtcslog)

        utime = runcmd(cmd=cmd, log=rtclog)
        return utime


    def run(self):
        label = self.label

        rushtime = self.rushver()
        oritime = self.oriver()

        resstr = "{}\tRealignerTargetCreator\trushtime:{}\toritime:{}".format(label, rushtime, oritime)
        return resstr


class IndelRealigner(object):
    def __init__(self, anadir, label):
        self.anadir = anadir
        self.label = label

    @staticmethod
    def __get_bam(bamdir):
        cancer, normal = "", ""
        for p, ds, fs in os.walk(bamdir):
            for f in fs:
                if f.endswith("_markdup.bam"):
                    if "case" in f or "cancer" in f:
                        cancer = os.path.join(p, f)
                    elif "normal" in f or "control" in f:
                        normal = os.path.join(p, f)
        if not cancer or not normal:
            logging.error("can't find bam in {}".format(bamdir))
            exit(1)
        return cancer, normal

    @staticmethod
    def __get_intervals(indir, label):
        respath = ""
        for p, ds, fs in os.walk(indir):
            for f in fs:
                if f == "coRealign.intervals":
                    if label in p:
                        respath = os.path.join(p, f)

        if not respath:
            logging.error("can't find RealignerTargetCreator results in {}".format(indir))
            exit(1)

        return respath

    def runshver(self):
        anadir = self.anadir
        bamdir = os.path.join(anadir, "output/picard")
        casebam, normalbam = IndelRealigner.__get_bam(bamdir)

        outdir = os.path.join(anadir, "output/realign/rush")
        logdir = os.path.join(anadir, "log/realign/rush")
        tempdir = os.path.join(anadir, "temp/realign/rush")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        rtcout = IndelRealigner.__get_intervals(indir=outdir, label="rush")
        logfile = os.path.join(logdir, "IndelRealigner.log")
        slog = os.path.join(logdir, "IndelRealigner.SUCCESS")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g  -Djava.io.tmpdir={temp} " \
              "-Djava.library.path={libpath} -jar {gatk} -T IndelRealigner -I {normalbam} -I {casebam} " \
              "--nWayOut _realign.bam -targetIntervals {rtcout} -known {mills} -known {dbsnp} -R {genome} " \
              "-allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment -dt NONE --maxReadsForRealignment 10000000 " \
              "-nThreads 32 > {log} 2>&1 && touch {slog} &&" \
              " mv ./*_sort_markdup_realign.b* {outdir}".format(java=java, temp=tempdir, libpath=gatklib,
                                                                normalbam=normalbam, casebam=casebam, rtcout=rtcout,
                                                                mills=mills, dbsnp=dbsnp, genome=genome, log=logfile,
                                                                slog=slog, gatk=newgatk, outdir=outdir)
        utime = runcmd(cmd=cmd, log=logfile)
        return utime

    def oriver(self):
        anadir = self.anadir
        bamdir = os.path.join(anadir, "output/picard")
        casebam, normalbam = IndelRealigner.__get_bam(bamdir)

        outdir = os.path.join(anadir, "output/realign/ori")
        logdir = os.path.join(anadir, "log/realign/ori")
        tempdir = os.path.join(anadir, "temp/realign/ori")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        rtcout = IndelRealigner.__get_intervals(indir=outdir, label="ori")
        logfile = os.path.join(logdir, "IndelRealigner.log")
        slog = os.path.join(logdir, "IndelRealigner.SUCCESS")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g  -Djava.io.tmpdir={temp} " \
              "-jar {gatk} -T IndelRealigner -I {normalbam} -I {casebam} " \
              "--nWayOut _realign.bam -targetIntervals {rtcout} -known {mills} -known {dbsnp} -R {genome} " \
              "-allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment -dt NONE --maxReadsForRealignment 10000000 " \
              "> {log} 2>&1 && touch {slog} &&" \
              " mv ./*_sort_markdup_realign.b* {outdir}".format(java=java, temp=tempdir, outdir=outdir,
                                                                normalbam=normalbam, casebam=casebam, rtcout=rtcout,
                                                                mills=mills, dbsnp=dbsnp, genome=genome, log=logfile,
                                                                slog=slog, gatk=oldgatk)
        utime = runcmd(cmd=cmd, log=logfile)
        return utime

    def run(self):
        label = self.label

        rushtime = self.runshver()
        oritime = self.oriver()

        resstr = "{}\tIndelRealigner\t{}\t{}".format(label, rushtime, oritime)

        return resstr


class BaseRecalibrator(object):
    def __init__(self, anadir, label):
        self.anadir = anadir
        self.label = label

    @staticmethod
    def __get_bam(bamdir):
        cancer, normal = "", ""
        for p, ds, fs in os.walk(bamdir):
            for f in fs:
                if f.endswith("_realign.bam"):
                    if "case" in f or "cancer" in f:
                        cancer = os.path.join(p, f)
                        cmd = "{} index {}".format(samtools, cancer)
                        if not os.path.isfile(cancer + ".bai"):
                            runcmd(cmd=cmd)
                    elif "normal" in f or "control" in f:
                        normal = os.path.join(p, f)
                        cmd = "{} index {}".format(samtools, normal)
                        if not os.path.isfile(normal + ".bai"):
                            runcmd(cmd=cmd)

        if not cancer or not normal:
            logging.error("can't find bam in {}".format(bamdir))
            exit(1)
        return cancer, normal

    def rushver(self, label):
        anadir = self.anadir
        bamdir = os.path.join(anadir, "output/realign/rush")
        casebam, normalbam = BaseRecalibrator.__get_bam(bamdir)

        outdir = os.path.join(anadir, "output/bqsr/rush")
        logdir = os.path.join(anadir, "log/bqsr/rush")
        tempdir = os.path.join(anadir, "temp/bqsr/rush")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        if label in ["case", "cancer"]:
            inbam = casebam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.case.log")
            slog = os.path.join(logdir, "BaseRecalibrator.case.SUCCESS")
        else:
            inbam = normalbam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.normal.log")
            slog = os.path.join(logdir, "BaseRecalibrator.normal.SUCCESS")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g " \
              "-Djava.io.tmpdir={temp} -Djava.library.path={lib} -jar {gatk} -T BaseRecalibrator " \
              "-I {inbam} -o {outgrp} -R {genome} -knownSites {mills} -knownSites {dbsnp} -knownSites {cosmic} -nct 32 " \
              "-allowPotentiallyMisencodedQuals -L {targetbed} -rbs 3000 > {log} 2>&1 && " \
              "touch {slog}".format(java=java, temp=tempdir, lib=gatklib, gatk=newgatk, inbam=inbam, outgrp=ougrp,
                                    genome=genome, mills=mills, dbsnp=dbsnp, cosmic=cosmic, log=logfile, slog=slog,
                                    targetbed=targetbed)

        utime = runcmd(cmd, logfile)
        return utime

    def oriver(self, label):
        anadir = self.anadir
        bamdir = os.path.join(anadir, "output/realign/ori")
        casebam, normalbam = BaseRecalibrator.__get_bam(bamdir)

        outdir = os.path.join(anadir, "output/bqsr/ori")
        logdir = os.path.join(anadir, "log/bqsr/ori")
        tempdir = os.path.join(anadir, "temp/bqsr/ori")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        if label in ["case", "cancer"]:
            inbam = casebam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.case.log")
            slog = os.path.join(logdir, "BaseRecalibrator.case.SUCCESS")
        else:
            inbam = normalbam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.normal.log")
            slog = os.path.join(logdir, "BaseRecalibrator.normal.SUCCESS")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g " \
              "-Djava.io.tmpdir={temp} -jar {gatk} -T BaseRecalibrator " \
              "-I {inbam} -o {outgrp} -R {genome} -knownSites {mills} -knownSites {dbsnp} -knownSites {cosmic} -nct 32 " \
              "-allowPotentiallyMisencodedQuals -L {targetbed} -rbs 3000 > {log} 2>&1 && " \
              "touch {slog}".format(java=java, temp=tempdir, gatk=oldgatk, inbam=inbam, outgrp=ougrp,
                                    genome=genome, mills=mills, dbsnp=dbsnp, cosmic=cosmic, log=logfile, slog=slog,
                                    targetbed=targetbed)

        utime = runcmd(cmd, logfile)
        return utime

    def run(self):
        label = self.label

        restr = []
        for type in ["case", "normal"]:
            rushtime = self.rushver(label=type)
            oritime = self.oriver(label=type)
            restr.append("{}_{}\tBaseRecalibrator\trushtime:{}\toritime:{}".format(label, type, rushtime, oritime))

        return "\n".join(restr)

class PrintReads(object):
    def __init__(self, anadir,label):
        self.anadir = anadir
        self.label = label

    @staticmethod
    def __get_bam(bamdir):
        cancer, normal = "", ""
        for p, ds, fs in os.walk(bamdir):
            for f in fs:
                if f.endswith("_realign.bam"):
                    if "case" in f or "cancer" in f:
                        cancer = os.path.join(p, f)
                        cmd = "{} index {}".format(samtools, cancer)
                        if not os.path.isfile(cancer + ".bai"):
                            runcmd(cmd=cmd)
                    elif "normal" in f or "control" in f:
                        normal = os.path.join(p, f)
                        cmd = "{} index {}".format(samtools, normal)
                        if not os.path.isfile(normal + ".bai"):
                            runcmd(cmd=cmd)

        if not cancer or not normal:
            logging.error("can't find bam in {}".format(bamdir))
            exit(1)
        return cancer, normal

    def rushver(self, label):
        anadir = self.anadir
        bamdir = os.path.join(anadir, "output/realign/rush")
        casebam, normalbam = PrintReads.__get_bam(bamdir)

        outdir = os.path.join(anadir, "output/bqsr/rush")
        logdir = os.path.join(anadir, "log/bqsr/rush")
        tempdir = os.path.join(anadir, "temp/bqsr/rush")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        if label in ["case", "cancer"]:
            inbam = casebam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.case.log")
            slog = os.path.join(logdir, "BaseRecalibrator.case.SUCCESS")
            outbam = os.path.join(outdir, "case_sort_markdup_realign_recal.bam")
        else:
            inbam = normalbam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.normal.log")
            slog = os.path.join(logdir, "BaseRecalibrator.normal.SUCCESS")
            outbam = os.path.join(outdir, "normal_sort_markdup_realign_recal.bam")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -XX:+UseParallelOldGC -Xms4g -Xmx16g " \
              "-Djava.io.tmpdir={temp} -Djava.library.path={lib} -jar {gatk} -T PrintReads -I {inbam} -BQSR {ougrp} " \
              "-o {outbam} -R {genome} -nct 32 -allowPotentiallyMisencodedQuals -rbs 3000 > {log} 2>&1 && " \
              "touch {slog} && {samtools} index {outbam}".format(java=java, temp=tempdir, lib=gatklib, gatk=newgatk,
                                                                 inbam=inbam, ougrp=ougrp, outbam=outbam, genome=genome,
                                                                 log=logfile, slog=slog, samtools=samtools)

        utime = runcmd(cmd=cmd, log=logfile)
        return utime

    def oriver(self, label):
        anadir = self.anadir
        bamdir = os.path.join(anadir, "output/realign/ori")
        casebam, normalbam = PrintReads.__get_bam(bamdir)

        outdir = os.path.join(anadir, "output/bqsr/ori")
        logdir = os.path.join(anadir, "log/bqsr/ori")
        tempdir = os.path.join(anadir, "temp/bqsr/ori")
        for subdir in [outdir, logdir, tempdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        if label in ["case", "cancer"]:
            inbam = casebam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.case.log")
            slog = os.path.join(logdir, "BaseRecalibrator.case.SUCCESS")
            outbam = os.path.join(outdir, "case_sort_markdup_realign_recal.bam")
        else:
            inbam = normalbam
            ougrp = os.path.join(outdir, "{}.grp".format(label))
            logfile = os.path.join(logdir, "BaseRecalibrator.normal.log")
            slog = os.path.join(logdir, "BaseRecalibrator.normal.SUCCESS")
            outbam = os.path.join(outdir, "normal_sort_markdup_realign_recal.bam")

        cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -XX:+UseParallelOldGC -Xms4g -Xmx16g " \
              "-Djava.io.tmpdir={temp} -jar {gatk} -T PrintReads -I {inbam} -BQSR {ougrp} " \
              "-o {outbam} -R {genome} -nct 32 -allowPotentiallyMisencodedQuals -rbs 3000 > {log} 2>&1 && " \
              "touch {slog} && {samtools} index {outbam}".format(java=java, temp=tempdir, gatk=oldgatk,
                                                                 inbam=inbam, ougrp=ougrp, outbam=outbam, genome=genome,
                                                                 log=logfile, slog=slog, samtools=samtools)

        utime = runcmd(cmd=cmd, log=logfile)
        return utime

    def run(self):
        label = self.label

        restr = []
        for type in ["case", "normal"]:
            rushtime = self.rushver(label=type)
            oritime = self.oriver(label=type)
            restr.append("{}_{}\tPrintReads\trushtime:{}\toritime:{}".format(label, type, rushtime, oritime))

        return "\n".join(restr)


def testgatk(anadir, label):
    s1 = RealignerTargetCreator(anadir=anadir, label=label).run()
    s2 = IndelRealigner(anadir=anadir, label=label).run()
    s3 = BaseRecalibrator(anadir=anadir, label=label).run()
    s4 = PrintReads(anadir=anadir, label=label).run()

    return "\n".join([s1, s2, s3, s4])

def opt():
    args = argparse.ArgumentParser(description="test gatk")
    args.add_argument("-a", "--anadir", help="The absolute path of anadir")
    args.add_argument("-l", "--label", help="The label of analysis")
    return args.parse_args()

def main():
    args = opt()
    anadir = args.anadir
    label = args.label
    setlog(logfile=os.path.join(anadir, "command.log"))

    ouf = os.path.join(anadir, "{}_time.log".format(label))

    with open(ouf, "w") as outfile:
        outfile.write(testgatk(anadir=anadir, label=label))

    logging.info("All analysis finished")

if __name__ == '__main__':
    main()
