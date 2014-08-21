#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Umstieg ROOT-Canvas -> matplotlib

Tag+Probe optional implementiert
"""

import matplotlib.pyplot as plt
import numpy as np

from scipy import optimize

from teilchen import Teilchen
import fitter

import math
import time

MIN_M_INV = 70
MAX_M_INV = 110

def particleFromFileLine(line):
	"""Get Muon data from formatted Line
	   Line-Format: <Pt> <Eta> <Phi> <M> <Charge> <Num Chambers> <Num Pixelhits> <Num Striphits> <Chi^2/n_DOF> <pfIso04> <sumPtIso03>"""
		
	line = line.split("\n")[0]
	args = line.split(" ")
	if len(args) == 4:
		# Neutrino
		return Teilchen(args[0], args[1], args[2], args[3])
	else:
		# Myon
		return Teilchen(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10])

def normalFill(diagramMassList, m1, m2):

	m = m1.invariantMass(m2)
	# volles Spektrum
	diagramMassList[9].append(m)
	
	# Gleiche Ladung? -> kein Z-Ereignis.
	if m1.charge() == m2.charge():
		diagramMassList[0].append(m)
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	elif m1.chi2nDOF() > 10 or m2.chi2nDOF() > 10:
		diagramMassList[2].append(m)

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	elif m1.numPixelHits() == 0 or m1.numStripHits() == 0 or \
		 m2.numPixelHits() == 0 or m2.numStripHits() == 0:
		diagramMassList[3].append(m)

	# Es müssen mindestens in 10 Kammern Hits sein
	elif m1.numChambers() <= 10 or m2.numChambers() <= 10:
		diagramMassList[4].append(m)
		
	# Rapidität kleiner 2,1
	elif m1.eta() > 2.1 or m2.eta() > 2.1:
		diagramMassList[5].append(m)
		
	# Gleicher Jet?
	elif m1.deltaR(m2) < 0.7:
		diagramMassList[1].append(m)

	# Mindestransversalimpuls 20 GeV
	elif m1.pt() < 20 or m2.pt() < 20:
		diagramMassList[6].append(m)

	# Ist eines der Myonen in einem Jet?
	elif m1.isolationFactor() > 1.15 or \
		 m2.isolationFactor() > 1.15:
		diagramMassList[7].append(m)

	elif m > MAX_M_INV or m < MIN_M_INV:
		diagramMassList[8].append(m)

	# Alle Tests bestanden!
	else:
		diagramMassList[10].append(m)

	return diagramMassList

def tagAndProbeFill(diagramMassList, m1, m2):

	m = m1.invariantMass(m2)

	diagramMassList[9].append(m) # volles Spektrum
	
	muon1valid = True
	muon2valid = True
	# Gleiche Ladung? -> kein Z-Ereignis.
	if m1.charge() == m2.charge():
		diagramMassList[0].append(m)
		return diagramMassList

	# Gleicher Jet?
	if m1.deltaR(m2) < 0.7:
		diagramMassList[1].append(m)
		return diagramMassList
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	if m1.chi2nDOF() > 10:
		muon1valid = False
	if m2.chi2nDOF() > 10:
		muon2valid = False

	if not muon1valid and muon2valid:
		diagramMassList[2].append(m)
		return diagramMassList

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	if m1.numPixelHits() == 0 or m1.numStripHits() == 0 and muon1valid:
		muon1valid = False
	if m2.numPixelHits() == 0 or m2.numStripHits() == 0 and muon2valid:
		muon2valid = False

	if not muon1valid and muon2valid:
		diagramMassList[3].append(m)
		return diagramMassList

	# Es müssen mindestens in 10 Kammern Hits sein
	if m1.numChambers() <= 10 and muon1valid:
		muon1valid = False
	if m2.numChambers() <= 10 and muon2valid:
		muon2valid = False

	if not muon1valid and muon2valid:
		diagramMassList[4].append(m)
		return diagramMassList
		
	# Pseudorapidität kleiner 2,1
	if m1.eta() > 2.1 and muon1valid:
		muon1valid = False
	if m2.eta() > 2.1 and muon2valid:
		muon2valid = False

	if not muon1valid and muon2valid:
		diagramMassList[5].append(m)
		return diagramMassList
	
	# Mindestransversalimpuls 20 GeV
	if m1.pt() < 20 and muon1valid:
		muon1valid = False
	if m2.pt() < 20 and muon2valid:
		muon2valid = False

	if not muon1valid and muon2valid:
		diagramMassListt[6].append(m)
		return diagramMassList

	if m1.isolationFactor() > 1.15 and muon1valid:
		muon1valid = False
	if m2.isolationFactor() > 1.15 and muon2valid:
		muon2valid = False

	if not (muon1valid and muon2valid):
		diagramMassList[7].append(m)
		return diagramMassList

	if m > MAX_M_INV or m < MIN_M_INV:
		diagramMassList[8].append(m)
		return diagramMassList

	# Mind. 1 Myon ist noch gültig -> Tag & Probe bestanden
	diagramMassList[10].append(m)

	return diagramMassList



def preParse(file, tagAndProbe):

	try:
		f = open(file,"r")
	except IOError:
		print "Datei wurde nicht gefunden."
		exit(0)

	#lineNums = getLen(file)/2 # 1 Event = 2 Zeilen
	lineNums = 4967428
	currentLine = 0
	currentPercent = 0
	startTime = int(time.time())
	lastNumber = 0
	lastSeconds = startTime
	print "Beginne. Parse %i Events"%lineNums

	l = f.readline()
	while l[0] == "#":
		l = f.readline()

	# DiagramMassList
	# Enthält alle Massen der gültigen Diagramme
	diagramMassList = []
	for i in range(11):
		diagramMassList.append([])
	
	while l:
		currentLine += 1
		if int(float(currentLine)*100/lineNums) != round(currentPercent):
			curRate = (currentLine - lastNumber)/(time.time()-lastSeconds)/1000
			avgRate = (currentLine)/(time.time()-startTime)/1000
			lastNumber = currentLine
			lastSeconds = time.time()
			currentPercent = round(float(currentLine)*100/lineNums)
			totalRate = currentLine/(time.time()-startTime)
			estimated = round((lineNums-currentLine)/totalRate)
			print "Habe %i Prozent geschafft. Current Rate: %.3f kHz, Avg: %.3f kHz, estimated remaining time: %is"%(currentPercent, curRate, avgRate, estimated)
			
		#if currentPercent == 5:
		#	break
		
		l1 = l
		l2 = f.readline()
		if not l2:
			raise IndexError("Inhalt der Eingabedatei ungültig. Zeilenanzahl muss Modulo3-Teilbar sein")

		m1 = particleFromFileLine(l1)
		m2 = particleFromFileLine(l2)

		if tagAndProbe: tagAndProbeFill(m1, m2)
		else: diagramMassList = normalFill(diagramMassList, m1, m2)

		# Nächste Zeile lesen
		l = f.readline()

	f.close()

	print "Parsing beendet. Benötigte Zeit: %i Sekunden. Avg Rate: %.3f kHz"%(int(time.time()-startTime), (currentLine)/(time.time()-startTime)/1000)
	drawPlots(diagramMassList)

	plt.figure(4)
	
	plt.title("Gefilterte Ereignisse")
	x = np.array([x/10.0 for x in range(MIN_M_INV*10, MAX_M_INV*10)])
	plt.hist(diagramMassList[10], x, color="b")
	plt.xlabel("$m_{inv}$")
	plt.ylabel("# Events")
	plt.xlim((MIN_M_INV, MAX_M_INV))

	f = open("diagrams-2/fits.txt", "w")
	# default gaussians
	p0 = np.array([-0.5, 1.0, 50.0, 2.0, 90.0])
	fitter.fitGauss(diagramMassList[10], f, p0,
					MAX_M_INV, MIN_M_INV, plt)
	f.write("\r\n")
	p0 = np.array([-0.05, 1.0, 90.0, 150.0, 2.0])
	fitter.fitBreitWiegner(diagramMassList[10], f, p0,
						   MAX_M_INV, MIN_M_INV, plt)
	f.close()
	plt.savefig("diagrams-2/zoomed.png")
	
	plt.show()

def drawPlots(diagramMassList, block = True):
	diagramsInfos = (("Ladungskriterium", "Myonen mit gleicher Ladung", 0, 0),
					 ("Richtungskriterium", "Myonen mit gleicher Richtung", 0, 1),
					 ("Spurqualität", "$\chi^2/n_{DOF}$ größer 10", 0, 2),
					 ("Detektorkriterium", "kein Signal im Spuren- oder Pixeldetektor", 1, 0),
					 ("Kammeranzahl", "Mindestens 10 Kammern", 1, 1),
					 ("Rapiditätskriterium", "Rapidität größer als 2,1", 1, 2),
					 ("Impulskriterium", "Impuls kleiner 20 GeV", 2, 0),
					 ("Isolationskriterium", "Myon ist nicht um Jet", 2, 1),
					 ("Massenfilter", "Invarianter Massenfilter", 2, 2),
					 ("Gesamtes Spektrum", "Gesamtes Myonen-Spektrum"),
					 ("Gefilterte Ereignisse", "Gefilterte Ereignisse"))

	d = 0.01
	f = 1.012
	binList = [0.3]
	for i in range(600):
		binList.append(binList[i]+d)
		d *= f

	f, ax = plt.subplots(3, 3)
	f.set_size_inches(12, 8, forward=True)
	f.subplots_adjust(left=0.125, right=0.9,
					  bottom=0.1, top=0.9,
					  wspace=0.4, hspace=0.8)

	for i in range(9):
		plotDiag(ax[diagramsInfos[i][2]][diagramsInfos[i][3]],
				 diagramMassList[i], binList,
				 diagramsInfos[i], True)

	plt.savefig("diagrams-2/filters.png")

	f = plt.figure(2)
	plotDiag(plt, diagramMassList[9], binList, diagramsInfos[9])
	plt.savefig("diagrams-2/spektrum.png")
	f = plt.figure(3)
	plotDiag(plt, diagramMassList[10], binList, diagramsInfos[10])
	plt.savefig("diagrams-2/filtered.png")

	#plt.show(block)

def plotDiag(to, data, binList, namesTup, isSubplot = False, log=True):
	if (log):
		to.loglog()
	if len(data) == 0 and isSubplot:
		to.set_title(namesTup[0])
		return
	if len(data) == 0 and not isSubplot:
		to.title(namesTup[0])
		return
	to.hist(data, binList, log = log, range=(0, 200), color="k", label="%i Events"%len(data))
	to.legend(frameon=False, fontsize='small')
	if isSubplot:
		to.set_xlabel("$m_{inv}$ [GeV]")
		to.set_ylabel("# Events")
		to.set_title(namesTup[0])
		to.set_xlim(1, 200)
	else:
		to.xlabel("$m_{inv}$ [GeV]")
		to.ylabel("# Events")
		to.title(namesTup[0])
		to.xlim(1, 200)

def getLen(file):
	f = open(file, "r")
	lineNum = 0
	l = f.readline()
	while l:
		l = f.readline()
		if l[0] == "#":
			continue
		lineNum += 1
	f.close()
	return lineNum

if __name__ == '__main__':
	f = "/mnt/usb/arbeit/2010-mit-iso/output.txt"
	preParse(f, False)
