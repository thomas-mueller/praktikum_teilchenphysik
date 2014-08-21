#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MuonFW5_Z.py
------------
Autor: Heinrich Blatt
Datum: 7.1.2013
Version: 1.0

Hauptprogramm zugehörig zu einer Bachelorarbeit.
Wertet die CMS-Myonendaten, die im Rahmen der
Bachelorarbeit für einen Versuch aufbereitet wurden,
aus. Generiert Diagramme für Z-Hosonen und fittet diese.

So könnte ein Beispielprogramm von Studenten für die
Auswertung von Z-Bosonen aussehen.

Dieses Programm ist sehr ausführlich, eine
minimale Version ist in der Datei
MuonFW5_Z-Minimal.py zu finden.
"""

# Definition der Teilchen-Klasse
from teilchen import Teilchen

# Die Module: fitter und plotter
import plotter

from fitpanel import Fitpanel

# Zur Geschwindigkeitsbestimmung und Fortschrittsdarstellung
import time

# Numpy um die Wurzel einer Liste / eines Arrays zu berechnen
import numpy as np

# Maximale und Minimale invariante masse
# Ist auch möglich direkt im Code, so ist er aber leichter
# auf W- und J/Psi-Mesonen wechselbar
MIN_M_INV = 83.5
MAX_M_INV = 97.5
NBINS = int((MAX_M_INV-MIN_M_INV)*1.5)

#eventsList = [open('../diagrams-2/events_%i.txt'%x, 'w') for x in range(12)]

def particleFromFileLine(line):
	"""Hole ein Myonenobjekt aus einer Datei

	   Line-Format: <Pt> <Theta> <Phi> <M> <Charge> \
<Num Chambers> <Num Pixelhits> <Num Striphits> \
<Chi^2/n_DOF> <pfIso04>"""
		
	line = line.split("\n")[0]
	args = line.split(",")
	if len(args) == 4:
		# Neutrino
		return Teilchen(args[0], args[1], args[2], args[3])
	else:
		# Myon
		return Teilchen(args[0], args[1], args[2], args[3], \
						args[4], args[5], args[6], args[7], \
						args[8], args[9])

def normalFill(fh, spektrum, dd, m1, m2):

	m = m1.invariantMass(m2)
	# volles Spektrum
	spektrum.fill(m)
	
	# Gleiche Ladung? -> kein Z-Ereignis.
	if m1.charge() == m2.charge():
		fh.fillSubdiagram('0', m)
		#eventsList[0].write(str(m1.evtPart())+",")
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	elif m1.chi2nDOF() > 10: # or m2.chi2nDOF() > 10:
		#eventsList[2].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('2', m)

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	elif m1.numPixelHits() == 0 or m1.numStripHits() == 0: # or \
		 #m2.numPixelHits() == 0 or m2.numStripHits() == 0:
		#eventsList[3].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('3', m)

	# Es müssen mindestens in 10 Kammern Hits sein
	elif m1.numChambers() <= 10: # or m2.numChambers() <= 10:
		#eventsList[4].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('4', m)
		
	# Rapidität kleiner 2,1
	
	#Die Triggereffizienz in den Endkappen (Bereiche > 2.1) bei dem genutzten
	#Algorithmus von 2010 nimmt stark ab. Es werden denoch Z's gefiltert, die
	#aber durch den Korrekturfaktor wieder berichtigt werden.
	#(CMS Collaboration, Performance of muon identification in 2010 data. 2011.
	# CMS PAS MUO-10-004)
	elif m1.eta() > 2.1: # or m2.eta() > 2.1:
		#eventsList[5].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('5', m)
		
	# Gleicher Jet?
	elif m1.deltaR(m2) < 0.7:
		#eventsList[1].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('1', m)

	# Ist eines der Myonen in einem Jet?
	elif m1.isolationFactor() > 1.15: # or \
		 #m2.isolationFactor() > 1.15:
		#eventsList[7].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('7', m)

	# Mindestransversalimpuls 20 GeV
	elif m1.pt() < 14: # or m2.pt() < 14:
		#eventsList[6].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('6', m)

	elif m > MAX_M_INV or m < MIN_M_INV:
		fh.fillSubdiagram('8', m)

	# Alle Tests bestanden!
	else:
		#eventsList[10].write(str(m1.evtPart())+",")
		dd.fill(m)

		#if m1.nTracks() == 2:
		#	eventsList[11].write(str(m1.evtPart())+",")

def tagAndProbeFill(fh, spektrum, dd, m1, m2):

	m = m1.invariantMass(m2)
	# volles Spektrum
	spektrum.fill(m)
	
	# Gleiche Ladung? -> kein Z-Ereignis.
	if m1.charge() == m2.charge():
		fh.fillSubdiagram('0', m)
		eventsList[0].write(str(m1.evtPart())+",")
		
	# Chi^2/nDOF-Wert (Güte des Ereignisses)
	elif m1.chi2nDOF() > 10:
		eventsList[2].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('2', m)

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	elif m1.numPixelHits() == 0 or m1.numStripHits() == 0:
		eventsList[3].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('3', m)

	# Es müssen mindestens in 10 Kammern Hits sein
	elif m1.numChambers() <= 10:
		eventsList[4].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('4', m)
		
	# Rapidität kleiner 2,1
	
	#Die Triggereffizienz in den Endkappen (Bereiche > 2.1) bei dem genutzten
	#Algorithmus von 2010 nimmt stark ab. Es werden denoch Z's gefiltert, die
	#aber durch den Korrekturfaktor wieder berichtigt werden.
	#(CMS Collaboration, Performance of muon identification in 2010 data. 2011.
	# CMS PAS MUO-10-004)
	#elif m1.eta() > 2.1:
	#	eventsList[5].write(str(m1.evtPart())+",")
	#	fh.fillSubdiagram('5', m)
		
	# Gleicher Jet?
	elif m1.deltaR(m2) < 0.7:
		eventsList[1].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('1', m)

	# Ist eines der Myonen in einem Jet?
	elif m1.isolationFactor() > 1.15:
		eventsList[7].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('7', m)

	# Mindestransversalimpuls 20 GeV
	elif m1.pt() < 14:
		eventsList[6].write(str(m1.evtPart())+",")
		fh.fillSubdiagram('6', m)

	elif m > MAX_M_INV or m < MIN_M_INV:
		fh.fillSubdiagram('8', m)

	# Alle Tests bestanden!
	else:
		eventsList[10].write(str(m1.evtPart())+",")
		dd.fill(m)

		if m1.nTracks() == 2:
			eventsList[11].write(str(m1.evtPart())+",")

def preParse(file, tagAndProbe):

	# Versuche die Datei zu öffnen
	try:
		f = open(file,"r")
	except IOError:
		print "Datei wurde nicht gefunden."
		exit(0)

	# optional: die Darstellung für den Fortschritt
	# (Initialisierung)
	lineNums = getLen(file)/2 # 1 Event = 2 Zeilen
	#lineNums = 4967428
	currentLine = 0
	currentPercent = 0
	startTime = int(time.time())
	lastNumber = 0
	lastSeconds = startTime
	print "Beginne. Parse %i Events"%lineNums


	fh = plotter.FilterHisto()
	diagrams = ((0, "Ladungskriterium", (0,0)), (1, "Richtungskriterium", (0,1)),
				(2, "Spurqualität", (0,2)),	 (3, "Detektorkritierum", (1,0)),
				(4, "Kammernzahl", (1,1)),	  (5, "Rapiditätskriterium", (1,2)),
				(6, "Impulskriterium", (2,0)),  (7, "Isolationskriterium", (2,1)),
				(8, "Massenfilter", (2,2)))
	for i in diagrams:
		fh.initSubdiagram(str(i[0]), i[1], i[2])

	spektrum = plotter.Histo("Spektrum")

	detaildiagram = plotter.DetailDiagram("Gefilterte Ereignisse",
										  MIN_M_INV, MAX_M_INV, NBINS)
	

	# Kopf der Datei wegspalten
	l = f.readline()
	while l[0] == "#":
		l = f.readline()

	# DiagramMassList
	# Generiere eine Liste mit 11 leeren Listen als Inhalt, die nachher
	# die invarianten Massen der gefilterten Variablen (die Filter und
	# auch die Filter passierenden) invarianten Massen.
	# In dieser Implementierung gibt es folgende Reihenfolge:
	# 1-9 Filter
	# 1. gleiche Ladung
	# 2. gleicher Jet
	# 3. Spurqualität
	# 4. Spurendetektor und Pixeldetektor-Hits vorhanden
	# 5. mindestens 10 Kammern
	# 6. Pseudorapidität kleiner 2,1
	# 7. Mindestransversalimpuls 20 GeV
	# 8. Isolationsfaktor des Myons
	# 9. Massenfilter (nur bei Tag+Probe)
	# 10. volles Spektrum
	# 11. gefiltertes Spektrum
	
	while l:
		# Fortschrittsanzeige - Optional
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
			
		if currentPercent == 1 and False:
			break

		# Zeile benennen und nächste Zeile auslesen
		l1 = l
		l2 = f.readline()

		# Generiert einen Fehler, falls die Dateiinhalte nicht Modulo2-Teilbar sind
		if not l2:
			raise IndexError("Inhalt der Eingabedatei ungültig. Zeilenanzahl muss Modulo3-Teilbar sein")

		# Erstelle Teilchen aus den Eingabezeilen
		m1 = particleFromFileLine(l1)
		m2 = particleFromFileLine(l2)

		# Filtere die Myonen nach den implementierten Filtern
		if tagAndProbe:
			tagAndProbeFill(fh, spektrum, detaildiagram, m1, m2)
		else:
			normalFill(fh, spektrum, detaildiagram, m1, m2)

		# Nächste Zeile lesen
		l = f.readline()

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
					MIN_M_INV, MAX_M_INV, NBINS, np.sqrt(binContents), '../diagrams-2/fits.txt')
	
	detaildiagram.plot()
	#detaildiagram.save("../diagrams-2/zoomed.png")

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
	f = "/portal/ekpcms5/home/tmueller/Praktikum_HBlatt/myhblatt/dimuon.txt"
	preParse(f, False)
