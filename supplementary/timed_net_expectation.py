# Analytic expectation of the timed Petri net cycle: delays T_i ~ N(mu_i, sd_i), the clamp, feed advance and
# feed return delays are multiplied by R = 0.5/P with P ~ N(0.5, sigma) truncated below at 0.28 and, for the
# completed-cycle statistics, conditioned on P >= 0.38.
import math
SEG = [(1.10,.11),(1.35,.14),(1.75,.15),(2.90,.24),(1.45,.17),(.90,.11),(.80,.09)]
SCALED = (0,3,4)
def moments(sigma, lo=0.38, n=200000):
    # numerical integration of E[R], E[R^2] for P ~ N(0.5, sigma) conditioned on P >= lo (lo >= 0.28 floor)
    a, b = lo, 0.5 + 8*sigma
    h = (b-a)/n
    def pdf(p): return math.exp(-0.5*((p-0.5)/sigma)**2)/(sigma*math.sqrt(2*math.pi))
    Z = sum(pdf(a+i*h)*(h if 0<i<n else h/2) for i in range(n+1))
    m1 = sum(pdf(a+i*h)*(0.5/(a+i*h))*(h if 0<i<n else h/2) for i in range(n+1))/Z
    m2 = sum(pdf(a+i*h)*(0.5/(a+i*h))**2*(h if 0<i<n else h/2) for i in range(n+1))/Z
    pcomplete = 0.5*(1+math.erf((0.5-lo)/(sigma*math.sqrt(2))))
    return m1, m2, pcomplete
for name, sigma in (('nominal',.028),('pressure variation',.075)):
    ER, ER2, pc = moments(sigma)
    muX = sum(SEG[i][0] for i in SCALED); varX = sum(SEG[i][1]**2 for i in SCALED)
    muY = sum(SEG[i][0] for i in range(7) if i not in SCALED); varY = sum(SEG[i][1]**2 for i in range(7) if i not in SCALED)
    mean = ER*muX + muY
    var = ER2*(varX+muX**2) - (ER*muX)**2 + varY
    print(name, 'E[R]=%.5f'%ER, 'P(complete)=%.4f'%pc, 'expected rejections of 300 = %.1f'%(300*(1-pc)),
          'mean=%.3f'%mean, 'sd=%.3f'%math.sqrt(var), 'unscaled sum', muY, 'scaled sum', muX)
print('sum of nominal means', sum(m for m,s in SEG), 'sd of plain sum', round(math.sqrt(sum(s*s for m,s in SEG)),3))
