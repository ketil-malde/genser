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

    mu = s/n
    var = ss/n-mu*mu
    return mu, sqrt(var)

from scipy.stats import norm

# assign a histogram to three distrs and re-estimate parameters
def expmax(distr, hist):
    h0 = {}
    h1 = {}
    h2 = {}    
    mu, sd, k0, k1, k2 = distr
    sd1, sd2, sd3 = sqrt(mu/2), sqrt(mu), sqrt(mu*2)
    # assign histogram to distrib
    for val,cnt in hist.items():
        pz = 0.0001  # uniform
        p0 = norm.pdf((val+0.5-mu/2)/sd1)  # prob of val under N(mu/2, sd/2)
        p1 = norm.pdf((val+0.5-mu  )/sd2)
        p2 = norm.pdf((val+0.5-mu*2)/sd3)
        ptot = p0 + p1 + p2 + pz
        h0[val] = cnt*p0/ptot
        h1[val] = cnt*p1/ptot
        h2[val] = cnt*p2/ptot
    # estimate the parameters
    mu0, sd0 = estimate(h0)
    mu1, sd1 = estimate(h1)
    mu2, sd2 = estimate(h2)
    print('  est:', (mu0,sd0), (mu1,sd1), (mu2,sd2))
    # weight estimate
    n0, n1, n2 = sum(h0.values()), sum(h1.values()), sum(h2.values())
    return ((mu0*2*n0 + mu1*n1 + mu2/2*n2)/(n0+n1+n2), (sd0*2*n0 + sd1*n1 + sd2*n2/2)/(n0+n1+n2), n0, n1, n2)

# criterion for end of convergence
def same(d0, d1):
    mu0, sd0, _, _, _ = d0
    mu1, sd1, _, _, _ = d1
    return(abs(mu0-mu1) < 0.005 and abs(sd1-sd0) < 0.005)

def errors(dist, hist):
    errs = {}
    mu, sd, k0, k1, k2 = dist
    for val, cnt in hist.items():
        e = cnt - (  k0 * norm.pdf( (val+0.5-mu/2)/(sd/2))
                   + k1 * norm.pdf( (val+0.5-mu)  / sd   )
                   + k2 * norm.pdf( (val+0.5-mu*2)/(sd*2)))
        errs[val] = e
    return errs

# iterate until convergence
# d0 = None
# while True:
#     d1 = expmax(d0, hist)
#    if same(d0,d1):
#        break
#    d0 = d1

# todo: write a test
#  - generate random data
#  - check that we can recover it
