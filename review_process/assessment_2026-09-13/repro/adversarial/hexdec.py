import sys, re
def hex2f(s):
    t = s.strip().replace(" ", "")
    neg = t.startswith("-")
    if t[0] in "+-": t = t[1:]
    mant, exp = t.split("X")
    v = float.fromhex("0x" + mant + "p" + str(int(exp, 16)))
    return -v if neg else v
for line in sys.stdin:
    print(re.sub(r'[+-][01]\.[0-9a-f]{13}X[+-][0-9a-f]{3}', lambda m: repr(hex2f(m.group(0))), line.rstrip()))
