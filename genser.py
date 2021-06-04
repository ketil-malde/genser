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
    if n<1:
        raise()
    mu = s/n
    var = ss/n-mu*mu
    return mu, var

from scipy.stats import poisson

def splithist(distr, hist):
    h0 = {}
    h1 = {}
    h2 = {}
    h3 = {}
    h4 = {}
    mu, sd0, sd1, sd2, sd3, sd4, k0, k1, k2, k3, k4 = distr
    # assign histogram to distrib
    pz = (k0+k1+k2+k3+k4)/5000
    for val,cnt in hist.items():
        p0 = k0 * poisson.pmf(int(mu/2), val)  # prob of val under N(mu/2, sd/2)
        p1 = k1 * poisson.pmf(int(mu),   val)
        p2 = k2 * poisson.pmf(int(mu*2), val)
        p3 = k3 * poisson.pmf(int(mu*3), val)
        p4 = k4 * poisson.pmf(int(mu*4), val)
        ptot = p0 + p1 + p2 + p3 + p4 + pz
        h0[val] = cnt*p0/ptot
        h1[val] = cnt*p1/ptot
        h2[val] = cnt*p2/ptot
        h3[val] = cnt*p3/ptot
        h4[val] = cnt*p4/ptot
    return h0, h1, h2, h3, h4

# assign a histogram to three distrs and re-estimate parameters
def expmax(distr, hist):
    h0, h1, h2, h3, h4 = splithist(distr, hist)
    # re-estimate the parameters
    mu0, var0 = estimate(h0)
    mu1, var1 = estimate(h1)
    mu2, var2 = estimate(h2)
    mu3, var3 = estimate(h3)
    mu4, var4 = estimate(h4)
    # weight estimate
    n0, n1, n2, n3, n4 = sum(h0.values()), sum(h1.values()), sum(h2.values()), sum(h3.values()), sum(h4.values())
    new_mu = (n0*mu0*2 + n1*mu1 + n2*mu2/2)/(n0+n1+n2) # mu1 # (mu0*2*n0 + mu1*n1 + mu2/2*n2 + mu3/3*n3)/(n0+n1+n2+n3)
    print(f'estimated distrs:  {mu0:.1f}±{var0:.1f} {mu1:.1f}±{var1:.1f} {mu2:.1f}±{var2:.1f} {mu3:.1f}±{var3:.1f} {mu4:.1f}±{var4:.1f} \r', end='')
    return (new_mu, var0, var1, var2, var3, var4, n0, n1, n2, n3, n4)

# criterion for end of convergence
def same(d0, d1):
    mu0, sd0, _, _, _, _, _, _, _, _, _ = d0
    mu1, sd1, _, _, _, _, _, _, _, _, _ = d1
    return(abs(mu0-mu1) < 0.1 and abs(sqrt(sd1)-sqrt(sd0)) < 0.1)

def errors(dist, hist):
    raise("don't use: it's poisson now")
    errs = {}
    mu, sd0, sd1, sd2, sd3, sd4, k0, k1, k2, k3, k4 = dist
    pz = (k0+k1+k2+k3+k4)/1000
    for val, cnt in hist.items():
        e = cnt - (  k0/sd0 * norm.pmf( (val+0.5-mu/2)/sd0)
                   + k1/sd1 * norm.pdf( (val+0.5-mu)  /sd1)
                   + k2/sd2 * norm.pdf( (val+0.5-mu*2)/sd2)
                   + k3/sd3 * norm.pdf( (val+0.5-mu*3)/sd3)
                   + k4/sd4 * norm.pdf( (val+0.5-mu*4)/sd4)
                   + pz)
        errs[val] = e
    return errs

def integr(hist):
    total = 0
    for val, count in hist.items():
        total = total + val*count
    return total

def integrate(distr, hist):
    h0, h1, _, _, _ = splithist(distr, hist)
    hs = {}
    for val, cnt in hist.items():
        hs[val] = cnt - h0[val] - h1[val]
    haploid = integr(h0)
    diploid = integr(h1)
    repeats = integr(hs)
    return haploid, diploid, repeats # NB! raw counts, divide by mu

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
