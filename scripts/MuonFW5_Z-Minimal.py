#!/usr/bin/env python
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

Dieses Programm ist die Minimalversion, um die Auswertung
durchzführen.
"""

# Definition der Teilchen-Klasse
from teilchen import Teilchen

# Die Module: fitter und plotter
import fitter
import plotter

# Numpy um die Wurzel einer Liste / eines Arrays zu berechnen
import numpy as np

# Maximale und Minimale invariante masse
# Ist auch möglich direkt im Code, so ist er aber leichter
# auf W- und J/Psi-Mesonen wechselbar
MIN_M_INV = 75
MAX_M_INV = 105

# Todo:
# - particleFromFileLine: Eta -> Theta
# - Plotfunktion weiter auslagern in eine Plotlib

def particleFromFileLine(line):
	"""Hole ein Myonenobjekt aus einer Datei

	   Line-Format: <Pt> <Eta> <Phi> <M> <Charge> \
<Num Chambers> <Num Pixelhits> <Num Striphits> \
<Chi^2/n_DOF> <pfIso04> <sumPtIso03>"""
		
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

	if not (muon1valid and muon2valid):
		diagramMassList[2].append(m)
		return diagramMassList

	# Spurendetektor und Pixeldetektor müssen jeweils Hits haben
	if m1.numPixelHits() == 0 or m1.numStripHits() == 0 and muon1valid:
		muon1valid = False
	if m2.numPixelHits() == 0 or m2.numStripHits() == 0 and muon2valid:
		muon2valid = False

	if not (muon1valid and muon2valid):
		diagramMassList[3].append(m)
		return diagramMassList

	# Es müssen mindestens in 10 Kammern Hits sein
	if m1.numChambers() <= 10 and muon1valid:
		muon1valid = False
	if m2.numChambers() <= 10 and muon2valid:
		muon2valid = False

	if not (muon1valid and muon2valid):
		diagramMassList[4].append(m)
		return diagramMassList
		
	# Pseudorapidität kleiner 2,1
	if m1.eta() > 2.1 and muon1valid:
		muon1valid = False
	if m2.eta() > 2.1 and muon2valid:
		muon2valid = False

	if not (muon1valid and muon2valid):
		diagramMassList[5].append(m)
		return diagramMassList
	
	# Mindestransversalimpuls 20 GeV
	if m1.pt() < 20 and muon1valid:
		muon1valid = False
	if m2.pt() < 20 and muon2valid:
		muon2valid = False

	if not (muon1valid and muon2valid):
		diagramMassList[6].append(m)
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

	# Versuche die Datei zu öffnen
	try:
		f = open(file,"r")
	except IOError:
		print "Datei wurde nicht gefunden."
		exit(0)

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
	diagramMassList = []
	for i in range(11):
		diagramMassList.append([])
	
	while l:
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
			diagramMassList = tagAndProbeFill(diagramMassList, m1, m2)
		else:
			diagramMassList = normalFill(diagramMassList, m1, m2)

		# Nächste Zeile lesen
		l = f.readline()

	f.close()

	# Zeichnet die Plots (Spektrum, Filter und selektierte Myonen)
	plotter.drawFilterPlots(diagramMassList, '../diagrams-2')

	detaildiagram = plotter.DetailDiagram("Gefilterterte Ereignisse",
								  diagramMassList[10], MIN_M_INV, MAX_M_INV, 1)

	binContents = detaildiagram.getBinContents()
	detaildiagram.drawErrors(binContents, np.sqrt(binContents))

	fp = Fitpannel(detaildiagram, diagramMassList[10],
					MIN_M_INV, MAX_M_INV, 20, np.sqrt(binContents))
	
	detaildiagram.save("../diagrams-2/zoomed.png")
	detaildiagram.show()

preParse("/home/heinrich/arbeit/2010/output-dimuon.txt", True)
