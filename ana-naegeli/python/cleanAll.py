#!/usr/bin/python

import os
import subprocess

CMD="lcgHome.py"
DEL="lcg-del -b -D srmv2 -l"
DIR="ntuples"
HOME="srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/naegelic"

p = subprocess.Popen(CMD + " " + DIR, shell=True, stdout=subprocess.PIPE, close_fds=True)
topdirs = p.communicate()[0]
topdirs = topdirs.split()

for dir in topdirs:
    d = os.path.basename(dir)
    p = subprocess.Popen(CMD + " " + DIR + "/" + d, shell=True, stdout=subprocess.PIPE, close_fds=True)
    rootfiles = p.communicate()[0]
    rootfiles = rootfiles.split()
    for root in rootfiles:
        r = os.path.basename(root)
        print DEL + " " + HOME + "/" + DIR + "/" + d + "/" + r
        os.system(DEL + " " + HOME + "/" + DIR + "/" + d + "/" + r)
