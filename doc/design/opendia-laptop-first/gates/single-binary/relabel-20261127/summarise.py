"""Summarise relabel outputs: per record and arm, old vs corrected pooled bounds and median N@FDP_c 1 %."""
import re, sys, glob, os
d = sys.argv[1]
pat = re.compile(r'^(\w+): pooled q<=0.01 lower ([\d.]+) \(90 % ([\d.]+)-([\d.]+)\), combined ([\d.]+) \(90 % ([\d.]+)-([\d.]+)\); E (\d+), N (\d+); median N@FDP_c 1 % (\d+) \((\d+)-(\d+)\)')
def read(path):
    out, r = {}, None
    for l in open(path):
        m = re.search(r'r = ([\d.]+)', l)
        if m and r is None: r = float(m[1])
        m = pat.match(l.strip())
        if m: out[m[1]] = m.groups()
    return r, out
print('| Record | Set | Arm | r old / corrected | E old / corrected | lower old / corrected | combined old / corrected | median N@FDP_c 1 % old / corrected |')
print('|---|---|---|---|---|---|---|---|')
for name in ['odp', 'sib', 'rtr', 'gbt', 'v7m6', 'blend']:
    for s in ['D', 'H']:
        o, c = f'{d}/{name}_{s}_old.txt', f'{d}/{name}_{s}_corrected.txt'
        if not os.path.exists(o): continue
        ro, po = read(o); rc, pc = read(c)
        for arm in po:
            a, b = po[arm], pc[arm]
            print(f'| {name.upper()} | {s} | {arm} | {ro:.4f} / {rc:.4f} | {a[7]} / {b[7]} | {float(a[1])*100:.2f} / **{float(b[1])*100:.2f} %** | '
                  f'{float(a[4])*100:.2f} / **{float(b[4])*100:.2f} %** | {int(a[9]):,} / **{int(b[9]):,}** |')
print()
print('| DIA-NN | r old / corrected | E old / corrected | lower old / corrected | combined old / corrected | N@FDP_c 1 % old / corrected |')
print('|---|---|---|---|---|---|')
for s in ['D', 'H']:
    vals = []
    for lab in ['old', 'corrected']:
        t = open(f'{d}/diann_{s}_{lab}.txt').read()
        r = float(re.search(r'r = ([\d.]+)', t)[1])
        m = re.search(r'DIA-NN: q<=0.01 N (\d+), E (\d+), lower ([\d.]+) \(90 % ([\d.]+)-([\d.]+)\), combined ([\d.]+) \(90 % ([\d.]+)-([\d.]+)\); N@FDP_c 1 % (\d+)', t)
        vals.append((r, m.groups()))
    (ro, a), (rc, b) = vals
    print(f'| {s} | {ro:.4f} / {rc:.4f} | {a[1]} / {b[1]} | {float(a[2])*100:.2f} / **{float(b[2])*100:.2f} %** | '
          f'{float(a[5])*100:.2f} / **{float(b[5])*100:.2f} %** | {int(a[8]):,} / **{int(b[8]):,}** |')
