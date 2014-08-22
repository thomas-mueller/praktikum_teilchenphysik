
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import fitter

import math

"""
plotter.py

Das Plotter-Modul.

Autor: Heinrich Blatt
Datum: Dezember 2013

Stellt die drei Klassen DetailDiagram, Histo und FilterHisto zur Verfügung.
Diese stellen Histogramme mithilfe von Matplotlib dar.
Es wird nicht die Matplotlib-Histogramm-Klasse genutzt, da diese Listen mit
Werten erwartet. Da diese Listen im Arbeitsspeicher liegen, ist diese Methode
sehr speicherintensiv. Die benutzte Methode speichert nur wie viele Inhalte
in den Bins existieren.
"""

class DetailDiagram(object):
	"""
	DetailDiagram

	Diese Klasse stellt ein DetailDiagramm bereit. Es enthält Fehlerbalken und
	ein nicht logarithmisch-skaliertes Histogramm.
	"""
	def __init__(self, title, min_m, max_m, n_bins, figure = 2):
		"""Generiert ein neues Histogramm.

		title - Überschrift des Histogramms
		min_m - minimale Masse
		max_m - maximale Masse
		n_bins - Anzahl der Bins"""
		# Wählt das Diagramm aus
		plt.figure(figure)
		# Benennt das Histogramm
		plt.title(title)
		# Generiert nullen in den Bins
		self._bin_contents = {}

		# Definiert die Klassenvariablen
		self._n_bins = n_bins
		self._min_m = min_m
		self._max_m = max_m
		self._figure = figure
		# Die Schrittweite pro Bin - wird oft gebraucht.
		self._step = float(max_m-min_m)/n_bins

		self._b = {}

	def fill(self, value, color='b'):
		"""Trägt den Wert value in das Histogramm ein"""
		# Welches Bin?
		try:
			bin_num = math.floor(float(value-self._min_m)/self._step)
		except TypeError:
			print "Versuche den Wert %s einzutragen, der jedoch nicht in den Grenzen des Histogramms liegt."%value
			return
		# Gibt es die Bins schon?
		if not color in self._bin_contents:
			self._bin_contents[color] = [0 for x in range(self._n_bins)]
		# Das bin inkrementieren
		try:
			self._bin_contents[color][bin_num] += 1
		except IndexError:
			print "Versuche den Wert %s einzutragen, der jedoch nicht in den Grenzen des Histogramms liegt."%value
			

	def getBinContents(self):
		"""Gibt die Inhalte der Bins zurück"""
		bin_contents = [0 for x in range(self._n_bins)]
		for i in self._bin_contents.values():
			for j in range(len(i)):
				bin_contents[j] += i[j]
		return bin_contents

	def drawErrors(self, bins, errors):
		"""Zeichnet die Fehler in die Bins"""
		# Diagramm auswählen
		plt.figure(self._figure)
		# die x-Werte des Bins herausfinden
		x = [(x*self._step)+self._min_m for x in range(self._n_bins)]
		# die Fehler in die Mitte der Bins schieben
		x = np.array(x)+(self._max_m-self._min_m)/(self._n_bins*2.0)
		# Und jetzt können wir einzeichnen.
		plt.errorbar(x, bins, yerr=errors, fmt='o', color='#AAAAAA')

	def save(self, file):
		"""Speichert das Diagramm in der Datei file"""
		# Diagramm auswählen
		plt.figure(self._figure)
		# und speichern
		plt.savefig(file)

	def plot(self, xlabel="m$_{\mu\mu}$ [GeV]"):
		"""Zeichnet das Diagramm"""
		if len(self._bin_contents) == 0:
			print "Es gibt keine Inhalte, die eingetragen werden könnten."
			return
		# Diagramm auswählen
		plt.figure(self._figure)
		# Einzelne Balken als Histogrammbalken zeichnen
		offset = [0 for x in range(self._n_bins)]
		for i in self._bin_contents.keys():
			for j in range(self._n_bins):
				if self._bin_contents[i] == 0: continue
				self._b[i] = plt.bar(j*self._step+self._min_m,
						self._bin_contents[i][j], self._step,
						offset[j], color=i)
				offset[j] += self._bin_contents[i][j]
				

		# Benennung und Grenzen der Achsen
		plt.xlabel(xlabel)
		plt.ylabel("# Events")
		plt.xlim((self._min_m, self._max_m))
		plt.ylim((0, max(offset)*1.2))

	def getplotobject(self):
		"""Gibt das plot-Objekt zurück"""
		# (Wichtig für Anbindung an den Fitpannel)
		return plt

	def addWLegend(self, labels):

		handles = []
		for i in labels.keys():
			handles.append(self._b[i])
		
		plt.legend(handles, labels.values())

class Histo(object):
	"""
	Klasse, die ein normales Histogramm mit doppelt logarithmischen Skalen
	beschreibt.

	Außerdem wird sie in die Klasse FilterHisto vererbt, das das 3x3-Fenster
	an Subplots angibt. Diese beinhalten die Filterplots.
	"""

	def __init__(self, title = '', figure = 1):
		""" Erzeuge Binliste. Benennt das Histogramm auch"""

		# Das Diagramm-auswählen und deren Zahl als Klassenvariable speichern
		plt.figure(figure)
		self._figure_num = figure

		# Generiere das Binning
		d = 0.01
		f = math.pow(1.012, 1.0/2)
		self.binList = [0.3]
		for i in range(900):
			self.binList.append(self.binList[i]+d)
			d *= f

		# Alle Bins auf 0 setzen
		self.bin_content = [0 for x in range(len(self.binList))]

		# Den Titel des Histogramms festlegen
		plt.title(title)

	def _findBin(self, value, left = 0, right = 0):
		"""Gibt das Bin mit dem Wert value zurück"""
		# Rekursiv definiert hat eine bessere Laufzeit als
		# eine for-schleife, die jedes Element prüft.
		if right == 0:
			# Initialisierung.
			right = len(self.binList)-1

		middleElem = int((right-left)/2+left)
		if right-left == 1:
			return left
		elif value > self.binList[middleElem]:
			return self._findBin(value, middleElem, right)
		elif value < self.binList[middleElem]:
			return self._findBin(value, left, middleElem)
		else:
			return middleElem

	def fill(self, value):
		"""Trägt den Wert value in das Histogramm ein"""
		binNum = self._findBin(value)
		self.bin_content[binNum] += 1

	def plot(self, xlabel = "m$_{\mu\mu}$ [GeV]"):
		"""Zeichnet das Histogramm"""
		# Diagramm auswäheln
		plt.figure(self._figure_num)
		# Achsen doppeltlogarithmisch setzen
		plt.loglog()
		# Histogrammbalken zeichnen
		self._drawHist(plt, self.bin_content)
		# Achsen beschriften
		plt.xlabel(xlabel)
		plt.ylabel("# Events")
		# X-Achsen-Bereich festlegen
		plt.xlim(0, 150)
		# Y-Achsen-Bereich festlegen
		plt.ylim(10, 500000)

	def save(self, file):
		"""Speichert das Histogramm in der Datei file"""
		# Diagramm auswählen
		plt.figure(1)
		# und speichern
		plt.savefig(file)

	def addLabels(self):
		"""Schreibt Label in das Histogramm"""
		plt.figure(1)
		plt.text(0.65, 17000, r"$\omega/\rho$")
		plt.text(0.95, 18000, r"$\phi$")
		plt.text(2.7, 170000, r"$J/\Psi$")
		plt.text(3.6, 9500, r"$\Psi'$")
		plt.text(9, 48000, r"$\Upsilon$")
		plt.text(85, 1500, "$Z$")

	def _drawHist(self, to, data):
		"""Zeichnet letztendlich das Histogramm

		schreibt jedes Bin als Balken in ein leeres Diagramm
		"""
		
		for i in range(len(self.binList)-1):
			if data[i] == 0: continue
			to.bar(self.binList[i], data[i],
				   self.binList[i+1]-self.binList[i], bottom=1, linewidth=0)

class FilterHisto(Histo):
	"""
	FilterHisto

	Generiert eine Übersicht die 9 Histogramme beinhaltet. Eignet sich gut,
	um die Werte der Filter einzutragen.
	"""
	def __init__(self):
		"""Initialisiert die Tabelle mit Histogrammen"""
		super(FilterHisto, self).__init__('', 2)
		self.diagrams = {}

	def initSubdiagram(self, histName, title, position):
		# Kompabilität - KillMe!
		self.initSubhisto(histName, title, position)

	def initSubhisto(self, histName, title, position):
		"""Initialisiert das Teilhistogramm mit dem Namen histName

		diagramName - Name des Histogramms
		name - Name der als Titel angezeigt wird
		position - tupel mit (x, y), (0,0) bis (2,2)"""

		# Es kann max. 9 Diagramme geben (3x3)
		if len(self.diagrams) > 9:
			print "Diagram kann nicht initialisiert werden."
			print "Es werden maximal 9 Filterdiagramme unterstützt."
			raise IndexError

		# Diagrammeigenschaften als Klassenvariablen speichern
		self.diagrams[histName] = (title, position[0], position[1],
								  [0 for x in range(len(self.binList))])

	def fillSubdiagram(self, histName, value):
		# Kompabilität - KillMe!
		self.fillSubhisto(histName, value)

	def fillSubhisto(self, histName, value):
		"""Trägt in das Teilhistogramm histName den Wert value ein"""
		binNum = self._findBin(value)
		self.diagrams[histName][3][binNum] += 1

	def plot(self, xlabels="m$_{\mu\mu}$ [GeV]"):
		"""Zeichnet das Histogramm"""
		# Diagramm in drei Teile teilen und formatieren
		self.f, ax = plt.subplots(3, 3)
		self.f.set_size_inches(12, 8, forward=True)
		self.f.subplots_adjust(left=0.125, right=0.9,
						  bottom=0.1, top=0.9,
						  wspace=0.4, hspace=0.8)

		# Alle Diagramme befüllen
		for i in self.diagrams.values():
			# Einträge zählen
			entries = 0
			for j in i[3]: entries += j
			print "Zeichne Plot %s mit %i Einträgen."%(i[0], entries)
			if entries == 0:  continue
			# In welches Diagramm schreiben wir?
			to = ax[i[1]][i[2]]
			# Doppellogarithmisch
			to.loglog()
			# Einzeichnen
			self._drawHist(to, i[3])
			to.set_xlabel(xlabels)
			to.set_ylabel("# Events")
			# Benennung des Diagramms
			to.set_title(i[0])
			# Grenzen der X-Achse
			to.set_xlim(1, 200)

	def save(self, file):
		"""Speichert die Histogrammsammlung in der Datei file"""
		# Figure auswählen
		plt.figure(2)
		# Diagramm speichern
		plt.savefig(file)

DetailHistogram = DetailDiagram
		
def show():
	# Wrapper um Matplotlib-Show()
	plt.show()
