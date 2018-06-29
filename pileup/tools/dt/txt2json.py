#!/usr/bin/python

import sys

if len(sys.argv)!=2:
    sys.exit("Error provide path as a command line argument")

path = sys.argv[1]

print(path+'lumis.txt')
f = open(path+'lumis.txt', 'r')
g = open(path+'lumis.json', 'w')

stuff = {}
for line in f:
    twins = line.split(" ", 1)
    tmpstuff = sorted(list(set(list(map(int,twins[1].rstrip().lstrip().rstrip(']').lstrip('[').split(","))))))
    shortstuff = []
    prev = tmpstuff[0]
    start = prev
    for nxt in tmpstuff[1:]:
        if nxt==prev+1:
            prev = nxt
        else:
            shortstuff.append("["+str(start)+","+str(prev)+"]")
            prev = nxt
            start = nxt
    shortstuff.append("["+str(start)+","+str(prev)+"]")
    stuff[int(twins[0])] = "["+", ".join(shortstuff)+"]"

nxtprnt = ""
cumulate = ""
for run in sorted(stuff.keys()):
    cumulate += nxtprnt
    cumulate += '"' + str(run) + '": '
    cumulate += stuff[run]
    nxtprnt = ", "
g.write("{" + cumulate + "}")

f.close()
g.close()
