#!/usr/bin/python

import sys

if len(sys.argv)!=2:
    sys.exit("Error provide path as a command line argument")

path = sys.argv[1]

f = open(path+'lumis.txt', 'r')
g = open(path+'lumis.json', 'w')

init = False
for line in f:
    twins = line.split(",", 1)
    if not init:
        init = True
        g.write('{')
    else:
        g.write(",\n ")
    g.write('"' + twins[0] + '": ' + twins[1].rstrip())
g.write("}")

f.close()
g.close()
