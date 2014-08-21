
# -*- coding: utf-8 -*-

import random
import time

random.seed(time.time())

TEMPLATE="/home/heinrich/arbeit/diagrams-3/events_%i.txt"
l = [TEMPLATE%x for x in range(12)]

o = open("/home/heinrich/arbeit/output-3.txt", "w")
for i in l:
	print "Suche aus Datei", i
	f = open(i, "r")
	line = f.read()
	if len(line) == 0:
		o.write("\n")
		continue
	f.close()
	li = line.split(",")
	for j in range(5):
		o.write(random.choice(li)+",")
	o.write("\n")
	o.flush()
o.close()



"""
import subprocess

f = open('output-1.txt', 'r')
o = open('output-files-1.txt', 'w')
content = f.read()
lines = content.split("\n")
for i in lines:
	if len(i) == 0:
		o.write("\n")
		continue
	sections = i.split(",")
	for j in sections:
		if len(j) == 0: continue
		event, run, lumi = j.split(" ")
		print "cmd 1:", 'python das_client.py \
--query="dataset dataset=/Mu/*2010*/AOD run=%s\"'%run
		o1 = subprocess.check_output('python das_client.py \
--query="dataset dataset=/Mu/*2010*/AOD run=%s"'%run)
		possible_datasets = []
		for k in o1.split("\n"):
			if k[0:3] == "/Mu": possible_datasets.append(k)
		for ds in possible_datasets:
			print "cmd2:", 'python das_client.py \
--query="file dataset=%s run=%s lumi=%s'%(ds, run, lumi)
			o2 = subprocess.check_output('python das_client.py \
--query="file dataset=%s run=%s lumi=%s'%(ds, run, lumi))
			file_output = o2.split("\n")
			br = False
			for l in file_output:
				if l[0:6] == "/store":
					o.write(l, ",")
					br = True
					break
			if br: break
		o.write("\n")
		o.flush()

"""

			
