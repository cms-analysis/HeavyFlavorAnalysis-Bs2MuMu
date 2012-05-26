#!/usr//bin/env python

import sys;

files = sys.argv[1:]

for s in files:
    f = open(s);
    cont = f.read();
    cont = cont.strip();
    print s + ":\t" + cont
    f.close()
