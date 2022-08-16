#!/usr/bin/env python3

# ARGS: 1. phylip concat 2. species tree topology 3. output rates 4. output trees

import sys
import subprocess

with open(sys.argv[1], "r") as f:
    subprocess.call("rm " + sys.argv[3], shell=True)
    subprocess.call("rm " + sys.argv[4], shell=True)
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        n = int(line.split()[0])
        with open("temp.phylip", "w") as f2:
            f2.write(line)
            for i in range(n):
                f2.write(f.readline())
        subprocess.call("iqtree -s temp.phylip -t " + sys.argv[2] + " -m GTR+G10+I -n 0 -wsr -redo", shell=True)
        subprocess.call("cat temp.phylip.rate >> " + sys.argv[3], shell=True)
        subprocess.call("cat temp.phylip.treefile >> " + sys.argv[4], shell=True)