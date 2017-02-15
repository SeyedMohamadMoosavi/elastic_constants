#! /usr/bin/env python
import sys

if len(sys.argv)>2:
    tol=float(sys.argv[2])
else:
    tol=3.0
print(tol)
f = open(sys.argv[1], 'r')
flines = f.readlines()
cangs=[]
for line in flines:
    if '_cell_angle_alpha' in line:
        cangs.append(float(line.split()[1]))
    if '_cell_angle_beta' in line:
        cangs.append(float(line.split()[1]))
    if '_cell_angle_gamma' in line:
        cangs.append(float(line.split()[1]))

cond=True
for ang in cangs:
    if not ((ang > 90.0-tol) and (ang < 90.0+tol)):
        cond=False
if len(flines)< 10:
    cond=False

if cond:
    f_ortho=open('orthogonals','a')
    f_ortho.write('%s\n'%sys.argv[1])
else:
    print(sys.argv[1])
