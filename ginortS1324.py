#!/usr/bin/python

import re
string = str(input())
low = []
up = []
num = []
odd = []
even = []
final = []

for x in string:
    if(x.isupper()):
        up.append(x)
    elif(x.islower()):
        low.append(x)
    elif(re.match("[0-9]",x)):
        num.append(int(x))
    else:
        continue

for a in sorted(low):
    final.append(a)

for b in sorted(up):
    final.append(b)


for c in num:
    if(c % 2 != 0):
        odd.append(c)
    else:
        even.append(c)

for o in sorted(odd):
    final.append(str(o))

for e in sorted(even):
    final.append(str(e))

for var in final:
    print("".join(var), end = "")
