# Literature to consider:
#   https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-128 Poisson hierchical model
# 


from random import randrange, gauss
from math import sqrt

import genser as G
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

hist = from_table()
hist.pop(0, None)  # delete zero counts (often inflated)

# estimate and compare        

mu_, var_ = G.estimate(hist)
dist = mu_, 1000.0, 2000.0, 1000.0, 1000.0, 1000.0

def get_modal():
    cur_cnt = 0
    for val, count in hist.items():
        if count > cur_cnt or val < mu_/2:
            cur_cnt = count
        else:
            return val

m = get_modal()
g = G.integr(hist)/mu_/1e6
print(f'Initial 1C estimates (using mu={mu_:.1f}, var={var_:.1f}): {g:.1f} Mbp; (using mode={m}): {g:.1f} Mbp')

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
mu, k0, k1, k2, k3, k4 = dist2

print(f'mus    {mu/2:3.1f} {mu:3.1f} {mu*2:3.1f} {mu*3:3.1f} {mu*4:.1f}')
print(f'counts {int(k0):10} {int(k1):10} {int(k2):10} {int(k3):10} {int(k4):10}')
print('')

from scipy.stats import poisson

def err_dist():
  for k in range(1,int(2*mu)):
    e0 = k0 * poisson.pmf(k, int(mu/2))  # prob of val under N(mu/2, sd/2)
    e1 = k1 * poisson.pmf(k, int(mu))
    e2 = k2 * poisson.pmf(k, int(mu*2))
    e3 = k3 * poisson.pmf(k, int(mu*3))
    e4 = k4 * poisson.pmf(k, int(mu*4))
    err = int(hist[k]-(e0+e1+e2+e3+e4))
    if err < 0:
        bar = '*'*(int(10*(-err)/hist[k]))
    else:
        bar = '-'*(int(10*err/hist[k]))
    print(f'{k:03} {int(hist[k]):10} pred: {int(e0):10} {int(e1):10} {int(e2):10} {int(e3):10} {int(e4):10} err: {err:10}  {bar}')

import matplotlib.pyplot as plt

def res_plot():
    limit = int(mu*4)
    plt.xlabel('Coverage')
    plt.ylabel('Count')
    #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
    xs = list(hist.keys())[:limit]
    ys = list(hist.values())[:limit]
    h0 = [k0*poisson.pmf(x, int(mu/2)) for x in xs]
    h1 = [k1*poisson.pmf(x, int(mu))   for x in xs]
    h2 = [k2*poisson.pmf(x, int(mu*2)) for x in xs]
    res = []
    for x in range(limit):
        res.append(ys[x]-h0[x]-h1[x]-h2[x])
    plt.plot(xs, ys, label='Coverage', linewidth=2)
    plt.plot(xs, h0)
    plt.plot(xs, h1)
    plt.plot(xs, h2)
    plt.plot(xs, res, ':', linewidth=2, label='Residuals')    
    plt.grid(True)

    plt.show()    

res_plot()
# err_dist()
    
hap, dip, rest = G.integrate(dist2, hist)
print()
print(f'Estimated total sequence, {(hap*2+dip+rest)/mu/1e6:f} Mbp,\n      haploid {hap*2/mu/1e6:f} diploid: {(dip+rest)/mu/1e6:f} (repeats: {rest/mu/1e6:f})')
print(f'Estimated total DNA (1C): {(hap+dip+rest)/mu/1e6:f} Mbp')
print()
