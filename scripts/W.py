#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

import time

from teilchen import Teilchen
from fitpanel import Fitpanel
import plotter

MIN_M_INV = 30
MAX_M_INV = 100
NBINS = int((MAX_M_INV-MIN_M_INV)*0.5)

#eventsList = [open('../diagrams-3/events_%i.txt'%x, 'w') for x in range(14)]

def particleFromFileLine(line):
	"""Get Muon data from formatted Line
	   Line-Format: <Muon1> or <Neutrino>
	   Muon-Format: <Pt> <Theta> <Phi> <M> <Charge> <Num Chambers> <Num Pixelhits> <Num Striphits> <Chi^2/n_DOF> <pfIso04> <nEvent> <nRun> <nLumi> <nVertices> <nTracks>
	   Neutrino-Format: <Pt> <Phi> <M> of missing Parts"""
	   
		
	line = line.split("\n")[0]
	args = line.split(",")

	if len(args) == 3:
		# Neutrino
		return Teilchen(args[0], 0, args[1], args[2])
	else:
		# Myon
		return Teilchen(args[0], args[1], args[2], args[3], args[4],
						args[5], args[6], args[7], args[8], args[9])

def neutrinoFill(fh, spektrum, dd, m, n, chargeList):

	mass = m.transverseInvariantMass(n)
	# volles Spektrum
	spektrum.fill(mass)
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	#elif m.chi2nDOF() > 10:
	if m.chi2nDOF() > 10:
		fh.fillSubdiagram('2', mass)
		#eventsList[2].write(str(m.evtPart())+",")

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	elif m.numPixelHits() == 0 or m.numStripHits() == 0:
		fh.fillSubdiagram('3', mass)
		#eventsList[3].write(str(m.evtPart())+",")

	# Es müssen mindestens in 10 Kammern Hits sein
	elif m.numChambers() < 10:
		fh.fillSubdiagram('4', mass)
		#eventsList[4].write(str(m.evtPart())+",")
		
	# Rapidität kleiner 2,1
	elif m.eta() > 2.1:
		fh.fillSubdiagram('5', mass)
		#eventsList[5].write(str(m.evtPart())+",")
		
	# Gleicher Jet?
	elif m.deltaR(n) < 0.7:
		fh.fillSubdiagram('1', mass)
		#eventsList[1].write(str(m.evtPart())+",")

	# Mindestransversalimpuls 20 GeV
	elif m.pT() < 20 or n.pT() < 20:
	#	#eventsList[6].write(str(m.evtPart())+",")
		fh.fillSubdiagram('6', mass)

	# Ist das Myon in einem Jet?
	elif m.isolationFactor() > 1.15:
		fh.fillSubdiagram('7', mass)
		#eventsList[7].write(str(m.evtPart())+",")

	#elif m.getEntfVertex() > 0.55:
	#	fh.fillSubdiagram('0', mass)

	elif mass > MAX_M_INV or mass < MIN_M_INV:
		#eventsList[10].write(str(m.evtPart())+",")
		fh.fillSubdiagram('8', mass)

	# Alle Tests bestanden!
	else:
		#dd.fill(mass)

		#dd2.fill(m.pT())
		#dd3.fill(n.pT())

		#if m.nVertices() == 1:
		#	eventsList[11].write(str(m.evtPart())+",")

		if m.charge() == 1:
			chargeList[0] += 1
			dd.fill(mass, 'b')
			#eventsList[12].write(str(m.evtPart())+",")
		else:
			chargeList[1] += 1
			dd.fill(mass, 'r')
			#eventsList[13].write(str(m.evtPart())+",")

	return chargeList
		
def preParse(file, tagAndProbe):

	try:
		f = open(file,"r")
	except IOError:
		print "Datei wurde nicht gefunden."
		exit(0)

	#lineNums = getLen(file)/2
	#lineNums = 31449331
	lineNums = 31827166
	currentLine = 0
	currentPercent = 0
	startTime = int(time.time())
	lastNumber = 0
	lastSeconds = startTime
	print "Beginne. Parse %i Events"%lineNums


	fh = plotter.FilterHisto()
	diagrams = ((0, "Anzahl Tracks", (0,0)), (1, "Richtungskriterium", (0,1)),
				(2, "Spurqualität", (0,2)),	 (3, "Detektorkritierum", (1,0)),
				(4, "Kammernzahl", (1,1)),	  (5, "Rapiditätskriterium", (1,2)),
				(6, "Impulskriterium", (2,0)),  (7, "Isolationskriterium", (2,1)),
				(8, "Massenfilter", (2,2)))
	for i in diagrams:
		fh.initSubdiagram(str(i[0]), i[1], i[2])

	spektrum = plotter.Histo("Spektrum")

	detaildiagram = plotter.DetailDiagram("Gefilterte Ereignisse",
										  MIN_M_INV, MAX_M_INV, NBINS)

	#dd2 = plotter.DetailDiagram("Myonen",
	#									  0, 80, 80, 4)

	#dd3 = plotter.DetailDiagram("Neutrinos",
	#									  0, 80, 80, 5)
	

	l = f.readline()
	while l[0] == "#":
		l = f.readline()

	chargeList = [0,0]
	
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
			
		#
		#if currentPercent == 1: # and False:
		#	break

		l1 = l
		l2 = f.readline()
		if not l2:
			raise IndexError("Inhalt der Eingabedatei ungültig. Zeilenanzahl muss Modulo2-Teilbar sein")

		m = particleFromFileLine(l1)
		n = particleFromFileLine(l2)

		chargeList = neutrinoFill(fh, spektrum, detaildiagram, m, n, chargeList)

		# Nächste Zeile lesen
		l = f.readline()

	f.close()

	print "Parsing beendet. Benötigte Zeit: %i Sekunden. Avg Rate: %.3f kHz"%(int(time.time()-startTime), (currentLine)/(time.time()-startTime)/1000)

	fh.plot("$m_{Transversal}$ [GeV]")

	binContents = detaildiagram.getBinContents()
	print "Anzahl Bins:", len(binContents)
	detaildiagram.drawErrors(binContents, np.sqrt(binContents))

	detaildiagram.plot("$m_{Transversal}$ [GeV]")
	#detaildiagram.save("../diagrams-3/zoomed.png")
	detaildiagram.addWLegend({'b':'$W^+$-Bosonen', 'r':'$W^-$-Bosonen'})

	spektrum.plot("$m_{Transversal}$ [GeV]")

	fp = Fitpanel(detaildiagram, binContents,
					MIN_M_INV, MAX_M_INV, NBINS,
					np.sqrt(binContents)) #, '../diagrams-3/fits.txt')
	
	print "Habe", chargeList[0], "positive und", chargeList[1], "negative Myonen gemessen."

	plotter.show()


def getLen(file):
	f = open(file, "r")
	lineNum = 0
	l = f.readline()
	while l:
		l = f.readline()
		if len(l) > 0 and l[0] == "#":
			continue
		lineNum += 1
	f.close()
	return lineNum

if __name__ == '__main__':
	f = "/mnt/usb/arbeit/output-w.txt"
	preParse(f, False)
