import os
import sys
import varibales
import logging
import subprocess
import time
import argparse
import multiprocessing
import setlog

def osmonitor(pid, res):
    monitor = varibales.monitor
    cmd = "{monitor} {pid} 2 > {res}".format(monitor=monitor, pid=pid, res=res)
    logging.debug("[monitor] {}".format(pid))
    os.system(cmd)

def opt():
    args = argparse.ArgumentParser(description="comparing gatk apps with and without FPGA")
    args.add_argument("-a", "--anadir", help="The current analysis directory", dest="anadir", required=True)
    return args.parse_args()

def rushver(case, control, anadir, targetbed):
    logging.info("Running rush gatk for {}".format(targetbed))
    usedtime = 0
    stime = time.time()
    outdir = os.path.join(anadir, "output/rush")
    logdir = os.path.join(anadir, "log/rush")
    temp = os.path.join(anadir, "temp/rush")
    for subdir in [outdir, logdir, temp]:
        if not os.path.isdir(subdir):
            os.makedirs(subdir)
            logging.debug("{} is not existed, making...".format(subdir))

    label = os.path.splitext(os.path.basename(targetbed))[0]
    outfile = os.path.join(outdir, os.path.splitext(os.path.basename(case))[0] + "{}.vcf".format(label))
    logfile = os.path.join(logdir, os.path.splitext(os.path.basename(case))[0] + "{}.rushgatk.log".format(label))

    mutect = varibales.mutect
    lib = varibales.mutectlib
    java = varibales.java
    genome = varibales.genome
    cosmic = varibales.cosmic
    dbsnp = varibales.dbsnp

    cmd = "{java} -d64 -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xms4g -Xmx16g -Djava.io.tmpdir={temp} " \
          "-Djava.library.path={lib} -jar {mutect} -T MuTect3 -R {genome} -I:tumor {case} -I:normal {control} " \
          "-L {targetbed} --dbsnp {dbsnp} --cosmic {cosmic} -contamination 0 --max_alt_alleles_in_normal_count 3 " \
          "--max_alt_alleles_in_normal_qscore_sum 40 --max_alt_allele_in_normal_fraction 0.02 -dt NONE -o {out} -nct 1 " \
          "--pair_hmm_implementation FPGA > {logfile} 2>&1 && " \
          "touch {logdir}/SUCCESS".format(java=java, temp=temp, lib=lib, mutect=mutect, genome=genome,
                                          case=case, control=control, targetbed=targetbed, dbsnp=dbsnp,
                                          cosmic=cosmic, out=outfile, logfile=logfile, logdir=logdir)

    cmdfile = os.path.join(logdir, "mutect.cmd")
    with open(cmdfile, "w") as ouf:
        ouf.write(cmd)

    logging.debug("[cmd] {}".format(cmd))

    subp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pid = subp.pid
    logging.debug("[pid] {}".format(pid))
    subp.wait()

    if subp.returncode != 0:
        logging.error("error reported for cmd: {}, detail info is {}".format(cmd, subp.stderr.read()))
        exit(1)
    else:
        etime = time.time()
        usedtime = (etime - stime)/60
        logging.info("Running rush gatk for {} finished.".format(targetbed))

    return "RushVer\t{targetbed}\t{utime}".format(targetbed=targetbed, utime=usedtime), outfile

def oriver(case, control, anadir, targetbed):
    logging.info("Running Ori gatk for {}".format(targetbed))
    usedtime = 0
    stime = time.time()
    outdir = os.path.join(anadir, "output/ori")
    logdir = os.path.join(anadir, "log/ori")
    temp = os.path.join(anadir, "temp/ori")
    for subdir in [outdir, logdir, temp]:
        if not os.path.isdir(subdir):
            os.makedirs(subdir)
            logging.debug("{} is not existed, making...".format(subdir))

    label = os.path.splitext(os.path.basename(targetbed))[0]

    outfile = os.path.join(outdir, os.path.splitext(os.path.basename(case))[0] + "{}.vcf".format(label))
    logfile = os.path.join(logdir, os.path.splitext(os.path.basename(case))[0] + "{}.oldgatk.log".format(label))

    oldgatk = varibales.oldgatk
    java = varibales.java
    genome = varibales.genome
    cosmic = varibales.cosmic
    dbsnp = varibales.dbsnp

    cmd = "{java} -d64 -server -XX:+UseParallelOldGC -XX:ParallelGCThreads=8 -Xms4g -Xmx16g -Djava.io.tmpdir={temp} " \
          "-jar {oldgatk} -T MuTect2 -R {genome} -I:tumor {case} -I:normal {control} -L {targetbed} --dbsnp {dbsnp} " \
          "--cosmic {cosmic} -contamination 0 --max_alt_alleles_in_normal_count 3 " \
          "--max_alt_alleles_in_normal_qscore_sum 40 --max_alt_allele_in_normal_fraction 0.02 -dt NONE " \
          "-o {outfile} > {logfile} 2>&1 && " \
          "touch {logdir}/SUCCESS".format(java=java, temp=temp, oldgatk=oldgatk, case=case, control=control,
                                          targetbed=targetbed, dbsnp=dbsnp, cosmic=cosmic, outfile=outfile,
                                          logfile=logfile, genome=genome, logdir=logdir)

    logging.debug("[cmd] {}".format(cmd))

    cmdfile = os.path.join(logdir, "mutect.cmd")
    with open(cmdfile, "w") as ouf:
        ouf.write(cmd)

    subp = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    pid = subp.pid
    logging.debug("[pid] {}".format(pid))
    subp.wait()

    if subp.returncode != 0:
        logging.error('error reported for cmd: {}, detail information is {}'.format(cmd, subp.stderr.read()))
        exit(1)
    else:
        etime = time.time()
        usedtime = (etime - stime)/60
        logging.info("Running Ori gatk for {} finished".format(targetbed))

    return "OriVer\t{targetbed}\t{utime}".format(targetbed=targetbed, utime=usedtime), outfile

def getbam(bamdir):
    case, control = "", ""
    for p, ds, fs in os.walk(bamdir):
        for f in fs:
            if f.endswith(".bam"):
                if "cancer" in f:
                    case = os.path.join(p, f)
                elif "normal" in f:
                    control = os.path.join(p, f)
    if not case or not control:
        logging.error("can't find bam file in {}".format(bamdir))
        exit(1)

    return case, control
