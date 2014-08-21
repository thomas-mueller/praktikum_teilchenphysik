#!/usr/bin/python
# -*- coding: utf-8 -*-

from PyQt4 import QtGui
from fitpannelgui import Ui_fitpannel

from scipy import special

import numpy as np
import ROOT as root

import math

"""
fitpannel.py

Das Fitpannel-Modul.

Autor: Heinrich Blatt
Datum: Januar 2014

Stellt die Fitpannel-Klasse zur Verfügung. Diese Funktioniert auf Basis von Qt.
Es ist daher sinnvoll, MatPlotLib's Qt4-Backend zu verwenden, da dann
man sich nicht der QApplication widmen muss.
Sollte man ein anderes Backend verwenden, so muss man die QApplication
initialisieren und starten, also einen eigenen Mainloop dafür generieren.
Dieser Mainloop muss in einem extra Thread gestartet werden, da er sonst
den Matplotlib-Mainloop blockiert.

Sollte man noch andere Funktionen anfügen wollen, so kann man diese leicht in
der Ui hinzufügen und hier als Funktion implementieren.
"""

gausstxt = """--- Gaußfit ---
Kurve gefittet. Parameter: a*x + b + Gauss.
a: %.5f err: %.5f, b: %.5f err: %.5f
#Prozess: %.2f err: %.2f, mean: %.5f err: %.5f, sigma: %.5f err: %.5f"""

breitwignertxt = """--- Breit-Wigner Fit ---
Kurve gefittet. Parameter: a*x + b + Breit-Wigner.
a: %.5f err: %.5f, b: %.5f err: %.5f
#Prozess: %.2f err: %.2f, mean: %.5f err: %.5f, width: %.5f err: %.5f"""

voigtiantxt = """--- Voigtianfit ---
Kurve gefittet. Parameter: a*x + b + Voigtian.
a: %.5f err: %.2f, b: %.5f err: %.5f
#Prozess: %.2f err: %.2f, mean: %.5f err: %.5f, width: %.5f err: %.5f, sigma: %.5f err: %.5f"""

class Fitpanel(QtGui.QMainWindow, Ui_fitpannel):
	"""
	Fitpannel-Klasse.

	Die eigentliche Gui-Generierung findet in der fitpannelgui.py statt, die
	aus der fitpannelgui.ui entstanden ist. Dieses ui-File wurde mit dem
	Qt4-Designer generiert. Mit diesem lässt es sich auch leicht editieren.
	"""

	def __init__(self, plotter, contents,
				 min_m, max_m, n_bins, yerrors, outputfile=''):
		"""
		Initialisiert den Fitpannel.

		plotter - plotter-Objekt um den Fit zu zeichnen
		contents - die Inhalte der Bins
		min_m, max_m, n_bins - Binning - Wichtig worauf gefittet werden soll
		yerrors - Fehler in y-Richtung
		outputfile - Optional, Datei in die das Ergebnis des Fits geschrieben wird
		"""
  
		# Parent-Initialisierung
		QtGui.QMainWindow.__init__(self)

		# Klassenvariablen
		self.plotter = plotter
		self._min_m = min_m
		self._max_m = max_m
		self._n_bins = n_bins
		self.plt = plotter.getplotobject()

		# Der letzte Fit wird in eine Datei geschrieben
		if outputfile != "": self.outputfile = open(outputfile, 'w')
		else: self.outputfile = ''

		# Zählt die Events die vorhanden sind, wichtig für Standardwert des
		# Signalstärke
		self.numEvents = 0
		for i in contents: self.numEvents += i

		# Plotobjekte - beinhaltet die Plots
		# [Gauss, Breitwigner, Voigtian]
		self.plotobjects = [False, False, False]
		self.fitresults = ['', '', '']

		# TH1F generieren aus dem die Fits generiert werden
		self.th = root.TH1F("invMass", "Invariant Mass", n_bins, min_m, max_m)
		self.th.SetFillColor(14)
		self.th.SetFillStyle(1001)
		
		for i in range(len(contents)):
			self.th.SetBinContent(i+1, contents[i])
			self.th.SetBinError(i+1, yerrors[i])
					
		# Das Userinterface initialisieren
		self.setupUi(self)

		# Verbinde die Buttons mit den Python-Funktionen
		self._connectButtons()

		# Fenster anzeigen
		self.show()

	def _connectButtons(self):
		"""
		Verbindet die Buttons mit den Funktionen
		"""
		self.do_gaussfit.clicked.connect(self._do_gaussfit)
		self.do_bwfit.clicked.connect(self._do_bwfit)
		self.do_voigtfit.clicked.connect(self._do_voigtfit)

		self.gauss_del_fit.clicked.connect(self._removegauss)
		self.bw_del_fit.clicked.connect(self._removebw)
		self.voigt_del_fit.clicked.connect(self._removevoigt)

	def __del__(self):
		"""
		Destruktor:
		Schreibt vor dem Löschen der Klasse noch die Fits in die Datei
		"""
		if self.outputfile != '':
			self.outputfile.write('\r\n'.join(self.fitresults))
			self.outputfile.close()

	def _removegauss(self):
		"""
		Löscht den Gaussfit
		"""
		if self.plotobjects[0] != False:
			self.plotobjects[0].set_linestyle('None')
			self.plt.draw()

	def _removebw(self):
		"""
		Löscht den Breit-Wigner-Fit
		"""
		if self.plotobjects[1] != False:
			self.plotobjects[1].set_linestyle('None')
			self.plt.draw()

	def _removevoigt(self):
		"""
		Löscht den Voigtianfit
		"""
		if self.plotobjects[2] != False:
			self.plotobjects[2].set_linestyle('None')
			self.plt.draw()
		
	def _dofit(self, pyFunc, rootfunc, funcname, write_labels, text,
			   plotnum, plotcolor, plotlabel, linestyle):

		"""
		Erstellt den eigentlichen Fit, schreibt die Ergebnisse in die
		Labels und zeichnet den Plot ein.
		"""
		#
		# 1.) Generiert den Fit mit Ergebnistext
		# 
		
		# Root erkennt Variablen anhand der internen Namen
		self.th.Fit(funcname, 'N+')

		# Alle Label sind mit Fehler (Anzahl: Labelzahl/2)
		# => Es gibt argCount argumente
		argCount = int(len(write_labels)/2)

		# Holt die Parameter und deren Fehler
		params = rootfunc.GetParameters()
		errors = rootfunc.GetParErrors()

		# ... Füllt eine Liste mit Parametern für den Ergebnistext und 
		# schreibt die Parameter und deren Fehler in den die das Widget
		textargs = []
		for i in range(argCount):
			textargs.append(float(params[i]))
			textargs.append(float(errors[i]))

			write_labels[i*2].setText("%.5f"%float(params[i]))
			write_labels[i*2+1].setText("%.5f"%float(errors[i]))

		# Gibt den Ergebnistext mit entsprechenden Werten aus und
		# speichert diesen in einer Kassenvariable, um ihn bei Beenden des
		# Fitpannels in eine Datei zu schreiben.
		print text%tuple(textargs)
		self.fitresults[plotnum-1] = text%tuple(textargs)			

		#
		# 2.) Zeichnet den Fit in das Diagramm ein.
		# 

		# Das Stepping ist 10 mal dichter als das Binning um glatte
		# Kurven zu gewährleisten
		stepping = float(self._max_m-self._min_m)/(self._n_bins*10)
		# generiere x-Werte
		x = np.array([x*stepping+self._min_m for x in range(self._n_bins*10)])
		# Unser Diagramm ist Diagramm 3 (Ggfs. Ändern!)
		self.plt.figure(2)
		
		if self.plotobjects[plotnum-1] == False:
			# Die Gaußfunktion ist noch nicht definiert
			# => Es wurde noch kein Plot eingezeichnet
			self.plotobjects[plotnum-1] = self.plt.plot(x, pyFunc(x, params, True),
										   color=plotcolor, label=plotlabel)[0]
			self.plt.legend(frameon=False, fontsize='small')
		else:
			# Es gibt schon einen Gaussplot
			# => Existierendem Plot neue Daten zuweisen
			self.plotobjects[plotnum-1].set_data(x, pyFunc(x, params, True))

		# Neuen Linienstil des Plots - Ggfs war der Plot schon gelöscht
		# (Linienstil = None)
		self.plotobjects[plotnum-1].set_linestyle(linestyle)

		# Und den Plot neu zeichnen ...
		self.plt.draw()

	def _do_gaussfit(self):
		"""Fittet gauss'sch aus den im Pannel angegebenen Daten."""

		# Definiere die Gaußfunktion als Python-Funktion
		step = float(self._max_m-self._min_m)/self._n_bins
		def gaussFunc(x, args, useNp = False):
			if not useNp: x = float(x[0])
			a, b = float(args[0]), float(args[1])
			sig, mean, sigma = float(args[2]), float(args[3]), float(args[4])
			if sig <= 0 or sigma <= 0: return 0
			return a*(x-mean) + b + \
				sig*step/(sigma*np.sqrt(2*np.pi))* \
				np.exp(-((x-mean)**2/(2*sigma**2)))

		# Importiere die Funktion in eine Root-Funktion
		self.gauss = root.TF1("gaussian", gaussFunc,
					   self._min_m, self._max_m, 5)

		# Parameter-(Label-)Liste generieren
		labels = (self.gauss_a, self.gauss_a_min, self.gauss_a_max,
				  self.gauss_b, self.gauss_b_min, self.gauss_b_max,
				  self.gauss_mean, self.gauss_mean_min, self.gauss_mean_max,
				  self.gauss_sigma, self.gauss_sigma_min, self.gauss_sigma_max)
		values = [float(x.text()) for x in labels]

		# Parametergrenzen setzen
		self.gauss.SetParameters(values[0], values[3],
							   0.8*self.numEvents, values[6], values[9])

		# Parameter setzen.
		# Bei Param 2 gibts ein Shift wegen dem Singalanteil
		# (Normierung des Gauss)
		for i in range(4):
			if i >= 2: val_i = i + 1
			else: val_i = i
			self.gauss.SetParLimits(val_i, values[i*3+1], values[i*3+2])
		
		# Ergebnislabels - werden von der _dofit-Funktion geschrieben
		write_labels = (self.gauss_a, self.gauss_a_err,
						self.gauss_b, self.gauss_b_err,
						self.gauss_num_signal, self.gauss_num_signal_err,
						self.gauss_mean, self.gauss_mean_err,
						self.gauss_sigma, self.gauss_sigma_err)

		# Jetzt: Übergabe von allem an die eigentliche Fitfunktion.
		self._dofit(gaussFunc, self.gauss, "gaussian",
					write_labels, gausstxt, 1, "g", "Gauß", 'solid')
		

	def _do_bwfit(self):
		"""Fittet breit-wiegner'sch aus den im Pannel angegebenen Daten."""

		# Definiere die Breit-Wigner-Funktion als Python-Funktion
		step = float(self._max_m-self._min_m)/self._n_bins
		def breitWigner(x, args, useNp = False):
			# params = (ax+b + breit wiegner)
			# breit wiegner = m, Max gamma
			if not useNp: x = float(x[0])
			a, b = float(args[0]), float(args[1])
			sig, m, width = float(args[2]), float(args[3]), float(args[4])
			if sig <= 0 or width <= 0: return 0
			return a*(x-m) + b + \
				   sig*step/((x-m)**2+0.25*width*width)

		# "Importiere" die Funktion in ROOT
		self.breitwigner = root.TF1("breitwigner", breitWigner,
									self._min_m, self._max_m, 5)

		# Parameter-(Label-)Liste generieren
		labels = (self.bw_a, self.bw_a_min, self.bw_a_max,
				  self.bw_b, self.bw_b_min, self.bw_b_max,
				  self.bw_mean, self.bw_mean_min, self.bw_mean_max,
				  self.bw_width, self.bw_width_min, self.bw_width_max)
		values = [float(x.text()) for x in labels]

		# Parameter-Standardwerte setzen
		self.breitwigner.SetParameters(values[0], values[3],
							   0.8*self.numEvents, values[6], values[9])

		# Parametergrenzen setzen.
		# Wieder Shift bei Parameter 2 wegen Signal
		for i in range(4):
			if i >= 2: val_i = i + 1
			else: val_i = i
			self.breitwigner.SetParLimits(val_i, values[i*3+1], values[i*3+2])


		# Ergebnislabels
		write_labels = (self.bw_a, self.bw_a_err,
						self.bw_b, self.bw_b_err,
						self.bw_num_signal, self.bw_num_signal_err,
						self.bw_mean, self.bw_mean_err,
						self.bw_width, self.bw_width_err)

		# ... Und der Fit.
		self._dofit(breitWigner, self.breitwigner, "breitwigner",
					write_labels, breitwignertxt, 2, "r", "Breit-Wigner", 'solid')
		

	def _do_voigtfit(self):
		"""Fittet die Voigtianfunktion aus den im Pannel angegebenen Daten."""

		# Definiere Python-Voigtianfunktion
		# Übersetzt in Python aus der RooFit-Implementation
		step = float(self._max_m-self._min_m)/self._n_bins
		voigtianConst = 1.0/np.sqrt(np.arctan2(0.0, -1.0))
		def voigtianFunc(x, args, useNp=False):
			if not useNp: x = float(x[0])
			a, b = float(args[0]), float(args[1])
			sig, mean = float(args[2]), float(args[3])
			w, s = float(args[4]), float(args[5])

			if sig <= 0 or s <= 0 or w <= 0:
				if useNp:   return [0 for x in range(len(x))]
				else:	   return 0

			coef = -0.5/(s*s)
			arg = x-mean

			# Voigtian
			c = 1.0/(np.sqrt(2.0)*s)
			a_ = 0.5*c*w
			u = c*arg

			# die RooFit implementation
			z = 0
			if useNp:   z = [complex(i, a_) for i in u]
			else:	   z = complex(u, a_)
			v = special.wofz(z)

			return a*(x-mean) + b + \
				   sig*step*c*voigtianConst*v.real

		# Und die Root-Voigitanfunktion
		self.voigtian = root.TF1("voigtian", voigtianFunc,
								 self._min_m, self._max_m, 6)

		# Parameter-(Label-)-Listen generieren
		labels = (self.voigt_a, self.voigt_a_min, self.voigt_a_max,
				  self.voigt_b, self.voigt_b_min, self.voigt_b_max,
				  self.voigt_mean, self.voigt_mean_min, self.voigt_mean_max,
				  self.voigt_width, self.voigt_width_min, self.voigt_width_max,
				  self.voigt_sigma, self.voigt_sigma_min, self.voigt_sigma_max)
		values = [float(x.text()) for x in labels]

		# Parameter-Standardwerte setzen
		self.voigtian.SetParameters(values[0], values[3], 0.8*self.numEvents,
									values[6], values[9], values[12])
									
		# Parametergrenzen setzen
		# Wieder ein Shift wegen Signalanteil
		for i in range(5):
			if i >= 2: val_i = i + 1
			else: val_i = i
			self.voigtian.SetParLimits(val_i, values[i*3+1], values[i*3+2])

		# Ergebnislabels
		write_labels = (self.voigt_a, self.voigt_a_err,
						self.voigt_b, self.voigt_b_err,
						self.voigt_num_signal, self.voigt_num_signal_err,
						self.voigt_mean, self.voigt_mean_err,
						self.voigt_width, self.voigt_width_err,
						self.voigt_sigma, self.voigt_sigma_err)

		# Und der Fit.
		self._dofit(voigtianFunc, self.voigtian, "voigtian",
					write_labels, voigtiantxt, 3, "y", "Voigt", 'solid')

