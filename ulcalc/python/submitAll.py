#!/usr/bin/python

import sys
import os

if len(sys.argv) < 3:
    print "usage: " + sys.argv[0] + " input.ulc prefix <additional>"
    exit(0)

INPUT = sys.argv[1]
OUTPUT = sys.argv[2]
DECAYS = ("bsmm","bdmm")
RESULTS = ("int","ul","sign")
BKG = "bkg"
ADD_ARGS = ""
if len(sys.argv) > 3:
    ADD_ARGS = sys.argv[3]

# Remove old stuff
print "==> Removing old directories..."
os.system("rm -rf " + OUTPUT)
os.system("mkdir " + OUTPUT)
os.chdir(OUTPUT)
os.system("mkdir " + BKG)
for res in RESULTS:
    os.system("mkdir %s" % res)

# Prepare the ulcs
print "==> Preparing ULC files..."
for j in range(0,100):
	os.system("cp ../%s %s/%s-%s-%d" % (INPUT,BKG,OUTPUT,BKG,j))
	for res in RESULTS:
		for dec in DECAYS:
			os.system("cp ../%s %s/%s-%s-%s-%d" % (INPUT,res,OUTPUT,dec,res,j))

print "==> Prepare Storage element..."
os.system("srmmkdir srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/naegelic/ulcalc/" + OUTPUT);

# Submit the jobs
BATCH_SCRIPT = os.environ['CMSSW_BASE'] + "/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/batch_calc.csh"
TAR = os.environ['CMSSW_BASE'] + "/src/HeavyFlavorAnalysis/Bs2MuMu/ulcalc/ulcalc.tar.gz"

os.chdir(BKG)
TOY_COUNT = 10000
ALGO_NAME = "bkg"
CMD = "run -q all.q -c %s -t %s -m batch -r 'STORAGE1 srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/naegelic/ulcalc/%s' -x 'bin/ulcalc --SM-exp --fixed-bkg --toys %d --seed %s -a %s' *%s*" % (BATCH_SCRIPT,TAR,OUTPUT,TOY_COUNT,ADD_ARG,ALGO_NAME,BKG)
print CMD
os.system(CMD)
os.chdir("..")

for res in RESULTS:
    for dec in DECAYS:
        # go to local directory
        os.chdir(res)
        # prepare special arguments to ulcalc depending on computation
        BDMM_OPT = ""
        if dec == "bdmm": BDMM_OPT = "--bdtomumu"
        TOY_COUNT = 1000;
        if res=="sign": TOY_COUNT = 10000
        if res=="sign": ALGO_NAME = "clb_hybrid"
        elif res=="ul": ALGO_NAME = "hybrid"
        else: ALGO_NAME = "int_hybrid"
        RANGE = 3
        if dec=="bdmm": RANGE *= 10
        RANGE_ARG = "-n 31 -r 0,%d" % RANGE
        if res=="sign": RANGE_ARG = ""
        CMD = "run -q all.q -c %s -t %s -m batch -r 'STORAGE1 srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/naegelic/ulcalc/%s' -x 'bin/ulcalc --SM-exp %s --fixed-bkg --toys %d --seed %s -a %s %s' *%s*" % (BATCH_SCRIPT,TAR,OUTPUT,BDMM_OPT,TOY_COUNT,ADD_ARGS,ALGO_NAME,RANGE_ARG,dec)
        print CMD
        os.system(CMD)
        os.chdir("..")
