# -*-coding:utf-8-*-

import os
import sys
import argparse
import logging

from getfqvol import FqVol
from setlog import setlog

def orivol(fqdir):
    vol = 0
    for p, ds, fs in os.walk(fqdir):
        for f in fs:
            if f.endswith("_1.fq.gz"):
                subvol = FqVol(fq=os.path.join(p, f)).value()
                vol += subvol
    return vol

class MagincDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def loadtimelog(timelog):
    resdict = {}
    with open(timelog, "r") as inf:
        for line in inf:
            fields = line.strip().split("\t")
            sample, anatype, rushtime, oritime = fields