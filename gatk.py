import os
import sys
import argparse
import logging
import subprocess

from setlog import setlog
from varibales import *

class RushVer(object):
    def __init__(self, case, control, anadir):
        self.case = case
        self.control = control
        self.anadir = anadir

    @staticmethod
    def __get_bam(anadir):
        casebam, controlbam = "", ""
        for p, ds, fs in os.walk(anadir):
            for f in fs:
                if f.endswith("_realign.bam"):
                    if "cancer" in f or "case" in f:
                        casebam = os.path.join(p, f)
                    elif "normal" in f or "control" in f:
                        controlbam = os.path.join(p, f)

        if not casebam or not controlbam:
            logging.error("can't find bam file in {}".format(anadir))
            exit(1)

        return casebam, controlbam

    @staticmethod
    def __local_realignment(case, control, anadir):
        caserealignbam, controlrealignbam = "", ""
        outdir = os.path.join(anadir, "output/realign")
        tempdir = os.path.join(anadir, "temp/realign")
        logdir = os.path.join(anadir, "log/realign")
        for subdir in [outdir, tempdir, logdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        rtclog = os.path.join(logdir, "RealignerTargetCreator.log")
        rtcslog = os.path.join(logdir, "RealignerTargetCreator.SUCCESS")
        illog = os.path.join(logdir, "IndelRealigner.log")
        ilslog = os.path.join(logdir, "IndelRealigner.SUCCESS")
        rtcout = os.path.join(outdir, "coRealign.intervals")

        rtccmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g  " \
              "-Djava.io.tmpdir={temp} -Djava.library.path={lib} -jar {gatk} -T RealignerTargetCreator " \
              "-I {normalbam} -I {casebam} -o {rtcout} -known {mills} -known {dbsnp} -L {targetbed} -R {genome} " \
              "-nt 12 -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment -dt NONE > {rtclog} 2>&1 && " \
              "touch {rtcslog}".format(java=java, temp=tempdir, lib=gatklib, gatk=newgatk, normalbam=control,
                                       casebam=case, rtcout=rtcout, mills=mills, dbsnp=dbsnp, targetbed=targetbed,
                                       genome=genome, rtclog=rtclog, rtcslog=rtcslog)

        ilcmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g  " \
                "-Djava.io.tmpdir={temp} -Djava.library.path={gatklib} -jar {gatk} -T IndelRealigner " \
                "-I {normalbam} -I {casebam} --nWayOut _realign.bam -targetIntervals {rtcout} -known {mills} " \
                "-known {dbsnp} -R {genome} -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment -dt NONE " \
                "--maxReadsForRealignment 10000000 -nThreads 32 > {illog} 2>&1 && " \
                "touch {ilslog} && " \
                "mv *_realign.* {outdir}".format(java=java, temp=tempdir, gatklib=gatklib, gatk=newgatk,
                                                 normalbam=control, casebam=case, rtcout=rtcout, mills=mills,
                                                 dbsnp=dbsnp, genome=genome,illog=illog, ilslog=ilslog, outdir=outdir)

        allcmd = "&&".join([rtccmd, ilcmd])
        logging.debug("[cmd] {}".format(allcmd))
        if os.path.isfile(ilslog):
            subp = subprocess.Popen(allcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subp.wait()

            if subp.returncode != 0:
                logging.error("Error reported for cmd: {}, detail info is {}".format(allcmd, subp.stderr.read()))
                exit(1)
            else:
                caserealignbam, controlrealignbam = RushVer.__get_bam(anadir)
        else:
            logging.info("{} existed, skipping local realignment...".format(ilslog))
            caserealignbam, controlrealignbam = RushVer.__get_bam(anadir)

        return caserealignbam, controlrealignbam


    @staticmethod
    def __bqsr(casebam, controlbam, anadir):

        outdir = os.path.join(anadir, "output/bqsr")
        tempdir = os.path.join(anadir, "temp/bqsr")
        logdir = os.path.join(anadir, "log/bqsr")
        for subdir in [outdir, tempdir, logdir]:
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

        brcmdstr = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g " \
                   "-Djava.io.tmpdir={tmp} -Djava.library.path={lib} -jar {gatk} -T BaseRecalibrator " \
                   "-I {bam} -o {grp} -R {genome} -knownSites {mills} -knownSites {dbsnp} -knownSites {cosmic} " \
                   "-nct 32 -allowPotentiallyMisencodedQuals -L {targetbed} -rbs 3000 > {brclog} 2>&1 && " \
                   "touch {brcslog}"

        prcmdstr = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -XX:+UseParallelOldGC -Xms4g -Xmx16g" \
                   " -Djava.io.tmpdir={tmp} -Djava.library.path={lib} -jar {gatk} -T PrintReads -I {bam} " \
                   "-BQSR {grp} -o {outbam} -R {genome} -nct 32 -allowPotentiallyMisencodedQuals -rbs 3000 > {log} 2>&1 && " \
                   "touch {slog}"

        casebrcgrp = os.path.join(outdir, "case.grp")
        casebrclog = os.path.join(logdir, "BaseRecalibrator.case.log")
        casebrcslog = os.path.join(logdir, "BaseRecalibrator.case.SUCCESS")
        caseprbam = os.path.join(outdir, "{}_recald.bam".format(os.path.splitext(os.path.basename(casebam))[0]))
        caseprlog = os.path.join(logdir, "PrintReads.case.log")
        caseprslog = os.path.join(logdir, "PrintReads.case.SUCCESS")
        casebrccmd = brcmdstr.format(java=java, tmp=tempdir, lib=gatklib, gatk=newgatk, bam=casebam, grp=casebrcgrp,
                                     genome=genome, mills=mills, dbsnp=dbsnp, cosmic=cosmic, targetbed=targetbed,
                                     brclog=casebrclog, brcslog=casebrcslog)
        caseprcmd = prcmdstr.format(java=java, tmp=tempdir, lib=gatklib, gatk=newgatk, bam=casebam, grp=casebrcgrp,
                                    outbam=caseprbam, genome=genome, log=caseprlog, slog=caseprslog)
        casecmd = "&&".join([casebrccmd, caseprcmd])

        controlbrcgrp = os.path.join(outdir, "control.grp")
        controlbrclog = os.path.join(logdir, "BaseRecalibrator.control.log")
        controlbrcslog = os.path.join(logdir, "BaseRecalibrator.control.SUCCESS")
        controlprbam = os.path.join(outdir, "{}_recald.bam".format(os.path.splitext(os.path.basename(controlbam))[0]))
        controlprlog = os.path.join(logdir, "PrintReads.control.log")
        controlprslog = os.path.join(logdir, "PrintReads.control.SUCCESS")
        controlbrccmd = brcmdstr.format(java=java, tmp=tempdir, lib=gatklib, gatk=newgatk, bam=controlbam,
                                        grp=controlbrcgrp, genome=genome, mills=mills, dbsnp=dbsnp, cosmic=cosmic,
                                        targetbed=targetbed, brclog=controlbrclog, brcslog=controlbrcslog)
        controlprcmd = prcmdstr.format(java=java, tmp=tempdir, lib=gatklib, gatk=newgatk, bam=controlbam, grp=controlbrcgrp,
                                       outbam=controlprbam, genome=genome, log=controlprlog, slog=controlprslog)
        controlcmd = "&&".join([controlbrccmd, controlprcmd])

        casesubp = subprocess.Popen(casecmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        controlsubp = subprocess.Popen(controlcmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        casesubp.wait()
        controlsubp.wait()

        logging.debug("[cmd] {}".format(casecmd))
        logging.debug("[cmd] {}".format(controlcmd))

        if casesubp.returncode != 0:
            logging.error("Error report for cmd: {}".format(casecmd))
            exit(1)

        if controlsubp.returncode != 0:
            logging.error("Error report for cmd: {}".format(controlcmd))
            exit(1)

    def run(self):
        case = self.case
        control = self.control
        anadir = self.anadir

        casebam, controlbam = self.__local_realignment(case=case, control=control, anadir=anadir)
        self.__bqsr(casebam=casebam, controlbam=controlbam, anadir=anadir)

def opt():
    args = argparse.ArgumentParser(description="run fpga gatk")
    args.add_argument("-a", "--anadir", help="The current analysis directory", dest="anadir", required=True)
    return args.parse_args()

def getbam(anadir):
    case, control = "", ""
    for p, ds, fs in os.walk(anadir):
        for f in fs:
            if f.endswith(".bam"):
                if "cancer" in f or "case" in f:
                    case = os.path.join(p, f)
                elif "normal" in f or "control" in f:
                    control = os.path.join(p, f)

    if not case or not control:
        logging.error("can't find bam file in {}".format(anadir))

    return case, control

def main():
    setlog()

    args = opt()
    anadir = args.anadir

    bamdir = os.path.join(anadir, "bam")

    case, control = getbam(bamdir)

    RushVer(case=case, control=control, anadir=anadir).run()

if __name__ == '__main__':
    main()

