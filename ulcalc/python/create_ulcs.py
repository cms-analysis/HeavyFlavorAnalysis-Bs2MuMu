#/usr/bin/python

import sys

if len(sys.argv) < 2 :
    print "usage: " + sys.argv[0] + " <template>"
    sys.exit(0)

#setup
template = sys.argv[1]
# tuples with first entry corresponding to barrel and
# second entry corresponding to endcap
bsmm_low = (0,0)
bsmm_high = (10,6)
bdmm_low = (0,0)
bdmm_high = (5,5)

print "Template file = " + template

# Open the template file
template_file = open(template)
template_str = template_file.read()
template_file.close()

for s0 in range(bsmm_low[0],bsmm_high[0]): # iterate through bsmm barrel
    for d0 in range(bdmm_low[0],bdmm_high[0]): # iterate through bdmm barrel
        for s1 in range(bsmm_low[1],bsmm_high[1]): #iterate through bsmm endcap
            for d1 in range(bdmm_low[1],bdmm_high[1]): #iterate through bdmm endcap
                t = template + "_bsb=" + str(s0) + "_bdb=" + str(d0) + "_bse=" + str(s1) + "_bde=" + str(d1)
                print "Saving '" + t + "'"
                f = open(t,"w")
                f.write(template_str)
                f.write("OBS_BSMM\t0\t" + str(s0) + "\n")
                f.write("OBS_BDMM\t0\t" + str(d0) + "\n")
                f.write("OBS_BSMM\t1\t" + str(s1) + "\n")
                f.write("OBS_BDMM\t1\t" + str(d1) + "\n")
                f.close()
