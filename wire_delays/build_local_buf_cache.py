"""Constructs the local buffer cache.
"""

import os
import sys
sys.path.insert(0, '..')

import tech

for node in tech.nodes:
    for N in [8]:
        os.system("python -u local_wires.py --K 6 --N %d --tech %s --density 0.5 --rebuffer 1" % (N, str(node)))
