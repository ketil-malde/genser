#!/usr/bin/python3
# Literature to consider:
#   https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-128 Poisson hierchical model
# 

# Parse args:
# -b bamfile, -t table
# -p --plot, -v --verbose

import argparse

p = argparse.ArgumentParser(description='Estimate genome size from coverage statistics.')
p.add_argument('-b', '--bamfile', action='store_true', help='Read coverage directly from BAM file')
p.add_argument('-p', '--plot', action='store_true', help='Generate coverage plot')
p.add_argument('-v', '--verbose', action='store_true', help='Show progress information')
p.add_argument('FILE', nargs='?', help='input file in tabular format (unless -b)')
args = p.parse_args()

# main program

from random import randrange, gauss
from math import sqrt

import util as G
import sys
import os
import subprocess

def run_command(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None)
    return iter(p.stdout.readline, b'')

def from_bamfile(filename):
    cmd = ['samtools', 'depth', '-d', '0', '-a', filename]
    hist = {}

    for l in run_command(cmd):
        i = int(l.split()[2])
        if i in hist:
            hist[i] += 1
        else:
            hist[i] = 1
    return hist

def from_table(fileobj):
    hist = {}
    for l in fileobj.readlines():
        ls = l.split()
        hist[int(ls[1])] = float(ls[0])
    return hist

if args.FILE is None:
    full_hist = from_table(sys.stdin)
elif args.bamfile:
    full_hist = from_bamfile(args.FILE)
else:
    with open(args.FILE, 'r') as fd:
        full_hist = from_table(fd)

# for efficiency, work on range 1..1000
hist = {}
for x in range(1,1000):
    if x in full_hist.keys():
        hist[x] = full_hist[x]

# estimate and compare        

mu_, var_ = G.estimate(hist)
r_, p_ = G.nbin_parms(mu_, var_)
dist = 2, r_, p_, 1000, 1000.0, 2000.0, 1000.0, 1000.0, 1000.0

def get_modal():
    cur_cnt = 0
    for val, count in hist.items():
        if count > cur_cnt or val < mu_/2:
            cur_cnt = count
        else:
            return val

m = get_modal()
g = G.integr(hist)/1e6
print(f'Initial 1C estimates (using mu={mu_:.1f}, var={var_:.1f}): {g/mu_:.1f} Mbp; (using mode={m}): {g/m:.1f} Mbp')

dist2 = G.expmax(dist, hist, args.verbose)
# print('**',dist2)

while(not G.same(dist,dist2)):
      dist = dist2
      dist2 = G.expmax(dist, hist, args.verbose)
#      es = G.errors(dist2, hist)
#      print('\r',dist2,end='')
#      print('  root sum square errs:', sqrt(sum([e*e for e in es.values()])))
#     print('  dist errs:', G.estimate(es))

mux, r, p, kx, k0, k1, k2, k3, k4 = dist2

if args.verbose:
    print(f'\n  counts (million) {int(kx)/1e6:.3f} {int(k0)/1e6:.3f} {int(k1)/1e6:.3f} {int(k2)/1e6:.3f} {int(k3)/1e6:.3f} {int(k4)/1e6:.3f}\n')

zs = full_hist[0]
full_hist.pop(0, None)
if args.plot:
    G.res_plot(full_hist, dist2)

low, hap, dip, rest = G.integrate(dist2, full_hist)
mu = r*p/(1-p)
print(f'Estimated total sequence, {(hap*2+dip+rest)/mu/1e6:f} Mbp,\n      haploid {hap*2/mu/1e6:f} diploid: {(dip+rest)/mu/1e6:f} (repeats: {rest/mu/1e6:f})')
print(f'Zero-coverage: {zs/1e6:.3f} Mbp, low coverage {low/1e6:.3f} Mbp')
print(f'Estimated total DNA (1C): {(hap+dip+rest)/mu/1e6:f} Mbp')

