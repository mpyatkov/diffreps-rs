"""Random inputs for comparing src/cdflib.rs with the DCDFLIB C code.

Lines: 'b a b x' -> cumbet(x, 1-x, a, b); 'n s r p' -> pnbinom(s, r, p).
"""
import math
import random

random.seed(7)


def lu(lo, hi):
    """Log-uniform number in [lo, hi]."""
    return math.exp(random.uniform(math.log(lo), math.log(hi)))


out = []
for _ in range(80000):  # generic incomplete beta
    a = lu(1e-8, 1e12)
    r = random.random()
    if r < .5:
        b = lu(1e-3, 1e5)
    elif r < .7:
        b = random.uniform(1, 2)
    elif r < .85:
        b = float(random.randint(1, 200))
    else:
        b = random.randint(1, 400) + 0.5
    t = random.random()
    x = random.random() if t < .5 else (1 - lu(1e-14, 1e-1) if t < .8 else lu(1e-14, 1e-1))
    out.append("b %r %r %r" % (a, b, x))

for _ in range(150000):  # negative binomial, as diffReps calls pnbinom
    m = lu(1, 5000)
    # variance at diffReps' floor (huge r, Poisson-like) or over-dispersed
    vm = m * m * 1e-8 * lu(1, 10) if random.random() < .3 else lu(1e-3, 100) * m
    r = m * m / vm
    p = r / (r + m)
    s = m + random.gauss(0, 3) * math.sqrt(m + vm)
    if random.random() < .3:
        s = math.floor(s) + random.choice([0.0, 0.5, 0.15106383, 0.84893617])
    if random.random() < .05:
        s = random.random()  # 0 < s < 1
    out.append("n %r %r %r" % (max(s, 0.0), r, p))

print("\n".join(out))
