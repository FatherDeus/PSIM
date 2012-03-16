#!/usr/bin/python

from sys import argv

file = open(argv[1],'r')
lines = list(file)
sortedlines = sorted(lines, key = lambda num: num[len(num)-2])
file.close()

file = open(argv[1],'w')
for s in sortedlines:
    file.write(s)
file.close()
