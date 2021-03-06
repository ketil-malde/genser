from math import sqrt

# mu is expected coverage for diploid, k0, k1, k2 is proportion (multipliers)
# for distributions around mu/2, mu, and 2*mu
# distr = None # (mu, k0, k1, k2)

# estimate a normal distribution for a histogram
def estimate(hist):
    n = s = ss = 0
    ls = sorted(hist.items())
    for cov, freq in ls[1:]:
        n += freq
        s += cov * freq
        ss += cov * cov * freq
    if n>=1:
        mu = s/n
        var = ss/n-mu*mu
        return mu, var
    else:
        return 0, 0

from scipy.stats import nbinom, poisson

def splithist(distr, hist, verbose=False):
    hx = {}
    h0 = {}
    h1 = {}
    h2 = {}
    h3 = {}
    h4 = {}
    mux, r, p, kx, k0, k1, k2, k3, k4 = distr
    # assign histogram to distrib
    pz = (kx+k0+k1+k2+k3+k4)/10000
    for val,cnt in hist.items():
        px = kx * poisson.pmf(val, mux)
        p0 = k0 * nbinom.pmf(val, r/2, 1-p)
        p1 = k1 * nbinom.pmf(val, r,   1-p)
        p2 = k2 * nbinom.pmf(val, r*2, 1-p)
        p3 = k3 * nbinom.pmf(val, r*3, 1-p)
        p4 = k4 * nbinom.pmf(val, r*4, 1-p)
        ptot = px + p0 + p1 + p2 + p3 + p4 + pz
        if ptot == pz:
            break
        hx[val] = cnt*px/ptot
        h0[val] = cnt*p0/ptot
        h1[val] = cnt*p1/ptot
        h2[val] = cnt*p2/ptot
        h3[val] = cnt*p3/ptot
        h4[val] = cnt*p4/ptot
#    if verbose:
#        print('*** dicts hX:', sum(h0.values()), sum(h1.values()), sum(h2.values()), sum(h3.values()))
    return hx, h0, h1, h2, h3, h4

def nbin_parms(mu, var):
    if var>mu*1.05:
        r = mu*mu/(var-mu)  # n is the number of successes (python)! r is failures (wikipedia)
        p = (var-mu)/var    # prob of any trial is a success
#        print(f'*** mu and var: {mu:.1f} {var:.1f}, r and p: {r:.2f}, {p:.2f}')
    else:
        r = mu*mu/(0.05*var)
        p = 0.05
    return (r, p)

# assign a histogram to three distrs and re-estimate parameters
def expmax(distr, hist, verbose=False):
    hx, h0, h1, h2, h3, h4 = splithist(distr, hist, verbose)
    # re-estimate the parameters
    mux, varx = estimate(hx)
    mu0, var0 = estimate(h0)
    mu1, var1 = estimate(h1)
    mu2, var2 = estimate(h2)
    mu3, var3 = estimate(h3)
    mu4, var4 = estimate(h4)
    # weight estimate
    nx, n0, n1, n2, n3, n4 = sum(hx.values()), sum(h0.values()), sum(h1.values()), sum(h2.values()), sum(h3.values()), sum(h4.values())
    mu = mu1 # (n0*mu0*2+n1*mu1+n2*mu2/2)/(n0+n1+n2)
    var = var1
    new_r, new_p = nbin_parms(mu, var)
    new_dist = ((mux+varx)/2, new_r, new_p, nx, n0, n1, n2, n3, n4)
    if verbose:
        print(f'  estimating distributions:  {mux:.1f}??{varx:.1f} {mu0:.1f}??{var0:.1f} {mu1:.1f}??{var1:.1f} {mu2:.1f}??{var2:.1f} {mu3:.1f}??{var3:.1f} {mu4:.1f}??{var4:.1f}\r', end='')
    return new_dist

# criterion for end of convergence
def same(d0, d1):
    mx0, r0, p0, _, _, _, _, _, _ = d0
    mx1, r1, p1, _, _, _, _, _, _ = d1
    return(abs(r0-r1) < 0.1 and abs(p0-p1) < 0.01)

def integr(hist):
    total = 0
    for val, count in hist.items():
        total = total + val*count
    return total

def integrate(distr, hist):
    hx, h0, h1, _, _, _ = splithist(distr, hist)
    hs = {}
    for val, cnt in hist.items():
        hs[val] = cnt - hx.get(val,0) - h0.get(val,0) - h1.get(val,0)
    lowcov  = integr(hx)
    haploid = integr(h0)
    diploid = integr(h1)
    repeats = integr(hs)
    return lowcov, haploid, diploid, repeats # NB! raw counts, divide by mu

import matplotlib.pyplot as plt

def res_plot(hist, distr):
    mux, r, p, kx, k0, k1, k2, _k3, _k4 = distr
    limit = 100 # int(mu*4)
    plt.xlabel('Coverage')
    plt.ylabel('Count')
    #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
    xs = list(hist.keys())[:limit]
    ys = list(hist.values())[:limit]
    hx = [kx*poisson.pmf(x, mux) for x in xs]
    h0 = [k0*nbinom.pmf(x, r/2, 1-p) for x in xs]
    h1 = [k1*nbinom.pmf(x, r,   1-p) for x in xs]
    h2 = [k2*nbinom.pmf(x, r*2, 1-p) for x in xs]
    res = []
    for x in range(limit):
        res.append(ys[x]-hx[x]-h0[x]-h1[x]-h2[x])
    plt.plot(xs, ys, '--', label='Coverage', linewidth=2)
    plt.plot(xs, hx)
    plt.plot(xs, h0)
    plt.plot(xs, h1)
    plt.plot(xs, h2)
    plt.plot(xs, res, ':', linewidth=2, label='Residuals')    
    plt.grid(True)

    plt.show()    

# Earlier version of res_plot, using textual output
def err_dist():
  for k in range(1,int(2*mu)):
    e0 = k0 * nbinom.pmf(k, r/2, 1-p)  # prob of val under N(mu/2, sd/2)
    e1 = k1 * nbinom.pmf(k, r,   1-p)
    e2 = k2 * nbinom.pmf(k, r*2, 1-p)
    e3 = k3 * nbinom.pmf(k, r*3, 1-p)
    e4 = k4 * nbinom.pmf(k, r*4, 1-p)
    err = int(hist[k]-(e0+e1+e2+e3+e4))
    if err < 0:
        bar = '*'*(int(10*(-err)/hist[k]))
    else:
        bar = '-'*(int(10*err/hist[k]))
    print(f'{k:03} {int(hist[k]):10} pred: {int(e0):10} {int(e1):10} {int(e2):10} {int(e3):10} {int(e4):10} err: {err:10}  {bar}')
