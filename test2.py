
from random import randrange, gauss
from math import sqrt

import genser as G
import sys

hist = {}
for l in sys.stdin.readlines():
    ls = l.split()
    hist[int(ls[1])] = float(ls[0])

# estimate and compare        

mu_, sd_ = G.estimate(hist) # wtf? tmp?
dist = mu_, sd_, sd_, sd_, sd_, sd_, 1000.0, 2000.0, 1000.0, 1000.0, 1000.0
# print(dist)

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
mu, sd0, sd1, sd2, sd3, sd4, k0, k1, k2, k3, k4 = dist2

print(f'mus    {mu/2:3.1f} {mu:3.1f} {mu*2:3.1f} {mu*3:3.1f} {mu*4:.1f}')
print(f'sigmas {sd0:3.1f} {sd1:3.1f} {sd2:3.1f} {sd3:3.1f} {sd4:3.1f}')
print(f'counts {int(k0):10} {int(k1):10} {int(k2):10} {int(k3):10} {int(k4):10}')
print('')

from scipy.stats import norm

for k in range(1,150):
    e0 = k0/sd0 * norm.pdf((k+0.5-mu/2)/sd0) # fuck: normalize by x width (area!)
    e1 = k1/sd1 * norm.pdf((k+0.5-mu)  /sd1)
    e2 = k2/sd2 * norm.pdf((k+0.5-mu*2)/sd2)
    e3 = k3/sd3 * norm.pdf((k+0.5-mu*3)/sd3)    
    e4 = k4/sd4 * norm.pdf((k+0.5-mu*3)/sd3)
    err = int(hist[k]-(e0+e1+e2+e3+e4))
    if err < 0:
        bar = '*'*(int(10*(-err)/hist[k]))
    else:
        bar = '-'*(int(10*err/hist[k]))
    print(f'{k:03} {int(hist[k]):10} pred: {int(e0):10} {int(e1):10} {int(e2):10} {int(e3):10} {int(e4):10} err: {err:10}  {bar}')
