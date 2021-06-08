# Literature to consider:
#   https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-128 Poisson hierchical model
# 


from random import randrange, gauss
from math import sqrt

import util as G
import sys
import subprocess

def run_command(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None)
    return iter(p.stdout.readline, b'')

def from_bamfile():
    cmd = ['samtools', 'depth', '-d', '0', '-a', sys.argv[1]]
    hist = {}

    for l in run_command(cmd):
        i = int(l.split()[2])
        if i in hist:
            hist[i] += 1
        else:
            hist[i] = 1
    return hist

def from_table():
    hist = {}
    for l in sys.stdin.readlines():
        ls = l.split()
        hist[int(ls[1])] = float(ls[0])
    return hist

full_hist = from_table()

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

dist2 = G.expmax(dist, hist)
# print('**',dist2)

while(not G.same(dist,dist2)):
      dist = dist2
      dist2 = G.expmax(dist, hist)
#      es = G.errors(dist2, hist)
#      print('\r',dist2,end='')
#      print('  root sum square errs:', sqrt(sum([e*e for e in es.values()])))
#     print('  dist errs:', G.estimate(es))

print('\n\nFinal:\n')
mux, r, p, kx, k0, k1, k2, k3, k4 = dist2

# print(f'mus    {mu/2:3.1f} {mu:3.1f} {mu*2:3.1f} {mu*3:3.1f} {mu*4:.1f}')
print(f'counts {int(kx):10} {int(k0):10} {int(k1):10} {int(k2):10} {int(k3):10} {int(k4):10}')
print('')

G.res_plot(full_hist, dist2)
# err_dist()

full_hist.pop(0,None)
hap, dip, rest = G.integrate(dist2, full_hist)
mu = r*p/(1-p)
print()
print(f'Estimated total sequence, {(hap*2+dip+rest)/mu/1e6:f} Mbp,\n      haploid {hap*2/mu/1e6:f} diploid: {(dip+rest)/mu/1e6:f} (repeats: {rest/mu/1e6:f})')
print(f'Estimated total DNA (1C): {(hap+dip+rest)/mu/1e6:f} Mbp')
print()
