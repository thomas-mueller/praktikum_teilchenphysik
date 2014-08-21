#!/usr/bin/python
# -*- coding: utf-8 -*-

""" 
MuonFW5_J_Psi.py
------------
Autor: Heinrich Blatt
Datum: 7.1.2013
Version: 1.0

Hauptprogramm zugehörig zu einer Bachelorarbeit.
Wertet die CMS-Myonendaten, die im Rahmen der
Bachelorarbeit für einen Versuch aufbereitet wurden,
aus. Generiert Diagramme für Z-Hosonen und fittet diese.

So könnte ein Beispielprogramm von Studenten für die
Auswertung von J/Psi-Mesonen aussehen.

Dieses Programm ist sehr ausführlich, eine
minimale Version ist in der Datei
MuonFW5_J_Psi-Minimal.py zu finden.
"""

# Definition der Teilchen-Klasse
from teilchen import Teilchen

# Die Module: fitter und plotter
from fitpanel import Fitpanel
import plotter

# Zur Geschwindigkeitsbestimmung und Fortschrittsdarstellung
import time

# Numpy um die Wurzel einer Liste / eines Arrays zu berechnen
import numpy as np

# Maximale und Minimale invariante masse
# Ist auch möglich direkt im Code, so ist er aber leichter
# auf W- und J/Psi-Mesonen wechselbar
#MIN_M_INV = 2.6
#MAX_M_INV = 3.6
MIN_M_INV=2.9
MAX_M_INV=3.3
NBINS = int((MAX_M_INV-MIN_M_INV)*10)*10

#eventsList = [open('../diagrams-1/events_%i.txt'%x, 'w') for x in range(12)]

def particleFromFileLine(line):
	"""Get Muon data from formatted Line
	   Line-Format: <Pt> <Theta> <Phi> <M> <Charge> \
<Num Chambers> <Num Pixelhits> <Num Striphits> \
<Chi^2/n_DOF> <pfIso04>"""
		
	line = line.split("\n")[0]
	args = line.split(",")
	
	return Teilchen(args[0], args[1], args[2], args[3],
					args[4], args[5], args[6], args[7],
					args[8], args[9])

def normalFill(fh, spektrum, dd, m1, m2):

	m = m1.invariantMass(m2)
	# volles Spektrum
	spektrum.fill(m)
	
	# Gleiche Ladung? -> kein J/Psi-Ereignis.
	if m1.charge() == m2.charge():
		fh.fillSubdiagram('0', m)
		#eventsList[0].write(str(m1.evtPart())+",")
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	elif m1.chi2nDOF() > 10: # or m2.chi2nDOF() > 10:
		fh.fillSubdiagram('2', m)
		#eventsList[2].write(str(m1.evtPart())+",")

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	elif m1.numPixelHits() == 0 or m1.numStripHits() == 0: # or \
		# m2.numPixelHits() == 0 or m2.numStripHits() == 0:
		fh.fillSubdiagram('3', m)
		#eventsList[3].write(str(m1.evtPart())+",")

	# Es müssen mindestens in 10 Kammern Hits sein
	elif m1.numChambers() <= 10: # or m2.numChambers() <= 10:
		fh.fillSubdiagram('4', m)
		#eventsList[4].write(str(m1.evtPart())+",")
		
	# Rapidität kleiner 2,1
	elif m1.eta() > 2.4: # or m2.eta() > 2.4:
		fh.fillSubdiagram('5', m)
		#eventsList[5].write(str(m1.evtPart())+",")
		
	# Gleicher Jet?
	elif m1.deltaR(m2) < 0.3:
		fh.fillSubdiagram('1', m)
		#eventsList[1].write(str(m1.evtPart())+",")

	# Mindestransversalimpuls 4 GeV
	elif m1.pt() + m2.pt() < 4 or m1.pt() + m2.pt() > 30 :
		fh.fillSubdiagram('6', m)
		#eventsList[6].write(str(m1.evtPart())+",")

	# Ist eines der Myonen in einem Jet?
	elif m1.isolationFactor() > 1.15: # or \
		 #m2.isolationFactor() > 1.15:
		fh.fillSubdiagram('7', m)
		#eventsList[7].write(str(m1.evtPart())+",")

	elif m > MAX_M_INV or m < MIN_M_INV:
		fh.fillSubdiagram('8', m)

	# Alle Tests bestanden!
	else:
		dd.fill(m)
		#eventsList[10].write(str(m1.evtPart())+",")

		#if m1.nTracks() == 2:
		#	eventsList[11].write(str(m1.evtPart())+",")

def tagAndProbeFill(diagramMassList, m1, m2):

	m = m1.invariantMass(m2)

	diagramMassList[9].append(m)
	
	muon1valid = True
	muon2valid = True
	# Gleiche Ladung? -> kein J/Psi-Ereignis.
	if m1.charge() == m2.charge():
		diagramMassList[0].append(m1)
		return diagramMassList

	# Gleicher Jet?
	if m1.deltaR(m2) < 0.7:
		diagramMassList[1].append(m1)
		return diagramMassList
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	if m1.chi2nDOF() > 10:
		muon1valid = False
	if m2.chi2nDOF() > 10:
		muon2valid = False

	if not (muon1valid or muon2valid):
		diagramMassList[2].append(m1)
		return diagramMassList

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	if m1.numPixelHits() == 0 or m1.numStripHits() == 0 and muon1valid:
		muon1valid = False
	if m2.numPixelHits() == 0 or m2.numStripHits() == 0 and muon2valid:
		muon2valid = False

	if not (muon1valid or muon2valid):
		diagramMassList[3].append(m1)
		return diagramMassList

	# Es müssen mindestens in 10 Kammern Hits sein
	if m1.numChambers() <= 10 and muon1valid:
		muon1valid = False
	if m2.numChambers() <= 10 and muon2valid:
		muon2valid = False

	if not (muon1valid or muon2valid):
		diagramMassList[4].append(m1)
		return diagramMassList
		
	# Pseudorapidität kleiner 2,1
	if m1.eta() > 2.1 and muon1valid:
		muon1valid = False
	if m2.eta() > 2.1 and muon2valid:
		muon2valid = False

	if not (muon1valid or muon2valid):
		diagramMassList[5].append(m1)
		return diagramMassList

	# Mindestransversalimpuls 20 GeV
	if m1.pt() < 20 and muon1valid:
		muon1valid = False
	if m2.pt() < 20 and muon2valid:
		muon2valid = False

	if not (muon1valid or muon2valid):
		diagramMassList[6].append(m1)
		return diagramMassList

	if m1.isolationFactor() > 1.15 and muon1valid:
		muon1valid = False
	if m2.isolationFactor() > 1.15 and muon2valid:
		muon2valid = False

	if not (muon1valid or muon2valid):
		diagramMassList[7].append(m)
		return diagramMassList

	if m > MAX_M_INV or m < MIN_M_INV:
		diagramMassList[8].append(m)
		return diagramMassList

	# Mind. 1 Myon ist noch gültig -> Tag & Probe bestanden
	diagramMassList[10].append(m1)
	return diagramMassList

def preParse(file, tagAndProbe):

	try:
		f = open(file,"r")
	except IOError:
		print "Datei wurde nicht gefunden."
		exit(0)

	#lineNums = getLen(file)/2 # 1 Event = 2 Zeilen
	#lineNums = 4967428
	lineNums = 4967427
	currentLine = 0
	currentPercent = 0
	startTime = int(time.time())
	lastNumber = 0
	lastSeconds = startTime
	print "Beginne. Parse %i Events"%lineNums

	l = f.readline()
	while l[0] == "#":
		l = f.readline()


	fh = plotter.FilterHisto()
	diagrams = ((0, "Ladungskriterium", (0,0)), (1, "Richtungskriterium", (0,1)),
				(2, "Spurqualität", (0,2)),	 (3, "Detektorkritierum", (1,0)),
				(4, "Kammernzahl", (1,1)),	  (5, "Rapiditätskriterium", (1,2)),
				(6, "Impulskriterium", (2,0)),  (7, "Isolationskriterium", (2,1)),
				(8, "Massenfilter", (2,2)))
	for i in diagrams:
		fh.initSubdiagram(str(i[0]), i[1], i[2])

	spektrum = plotter.Histo("Spektrum")

	detaildiagram = plotter.DetailDiagram("Gefilterterte Ereignisse",
										  MIN_M_INV, MAX_M_INV, NBINS)
	
	
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
			
		#if currentPercent == 10:
		#	break
		
		l1 = l
		l2 = f.readline()
		if not l2:
			raise IndexError("Inhalt der Eingabedatei ungültig. Zeilenanzahl muss Modulo3-Teilbar sein")
		

		m1 = particleFromFileLine(l1)
		m2 = particleFromFileLine(l2)

		if tagAndProbe: tagAndProbeFill(fh, spektrum, detaildiagram, m1, m2)
		else: normalFill(fh, spektrum, detaildiagram, m1, m2)

		# Nächste Zeile lesen
		l = f.readline()

	f.close()

	#for f in eventsList:
	#	f.close()

	# Statistik
	print "Parsing beendet. Benötigte Zeit: %i Sekunden. Avg Rate: %.3f kHz"%(int(time.time()-startTime), (currentLine)/(time.time()-startTime)/1000)
	print "Zeichne Plots ..."
	t = time.time()

	binContents = detaildiagram.getBinContents()
	detaildiagram.drawErrors(binContents, np.sqrt(binContents))

	print "Zeichnen beendet. Benötigte Zeit: %i Sekunden"%int(time.time()-t)

	fh.plot()

	spektrum.addLabels()
	spektrum.plot()

	fp = Fitpanel(detaildiagram, binContents,
					MIN_M_INV, MAX_M_INV, NBINS, np.sqrt(binContents), '../diagrams-1/fits.txt')
	
	#detaildiagram.save("../diagrams-1/zoomed.png")
	detaildiagram.plot()

	plotter.show()

def getLen(file):
	f = open(file, "r")
	lineNum = 0
	l = f.readline()
	while l:
		l = f.readline()
		if len(l) == 0 or l[0] == "#":
			continue
		lineNum += 1
	f.close()
	return lineNum

if __name__ == '__main__':
	f = "/mnt/usb/arbeit/output-dimuon.txt"
	preParse(f, False)
