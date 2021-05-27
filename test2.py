
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
dist = mu_, sd_, 1000.0, 1000.0, 1000.0, 1000.0
print(dist)

dist2 = G.expmax(dist, hist)
print('**',dist2)

while(not G.same(dist,dist2)):
      dist = dist2
      dist2 = G.expmax(dist, hist)
      print('\r',dist2,end='')
      es = G.errors(dist2, hist)
      print('  root sum square errs:', sqrt(sum([e*e for e in es.values()])))
#     print('  dist errs:', G.estimate(es))

print('\nFinal:')      
print(dist2)      

mu, sd, k0, k1, k2, k3 = dist

from scipy.stats import norm

for k in range(1,100):
    e1 = k0*norm.pdf((k+0.5-mu/2)/(sd/2))
    e2 = k1*norm.pdf((k+0.5-mu)/sd)
    e3 = k2*norm.pdf((k+0.5-mu*2)/(sd*2))
    print(k, hist[k], 'pred:', e1, e2, e3, 'err:', hist[k]-(e1+e2+e3))
