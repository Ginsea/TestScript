# -*- coding:utf-8 -*-

import os
import sys
import logging

cosmic = "/data/chenzhx/FPGA/db/known_db/COSMIC/b37_cosmic_v73_061615.vcf.gz"
targetbed = "/data/chenzhx/FPGA/config/nc_v13/regs/flk50_chip.bed"
dbsnp = "/data/chenzhx/FPGA/db/known_db/dbSNP/dbsnp_138.b37.del100.vcf.gz"
mills = "/data/chenzhx/FPGA/db/known_db/1000G_gold_standard/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
genome = "/data/chenzhx/FPGA/db/genome/hs37d5.fa"

mutectlib = "/data/chenzhx/FPGA/apps/FPGA/JNILib/"
mutect = "/data/chenzhx/FPGA/apps/FPGA/fpga.jar"
oldgatk = "/data/chenzhx/FPGA/apps/GenomeAnalysisTK.Ori.jar"
newgatk = "/data/chenzhx/FPGA/apps/GATK/gatk.jar"
gatklib = "/data/chenzhx/FPGA/apps/GATK/JNILib/"
java = "/usr/lib/jvm/java-1.8.0/bin/java"
bwa = "/data/chenzhx/FPGA/apps/bwaMod"
ncfilter="/data/chenzhx/FPGA/apps/NCfilter"
picard="/data/chenzhx/FPGA/apps/picard.jar"
samtools='/data/chenzhx/FPGA/apps/samtools'

monitor="/data/chenzhx/FPGA/apps/psMonitor.sh"