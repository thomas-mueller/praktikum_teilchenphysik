
# -*- coding: utf-8 -*-

import numpy as np
from scipy import optimize

from scipy import special

from ROOT import *

import math

# Startwerte der Fitvariablen
# Aufgebaut nach dem Schema
# <variable1>, <variable2>, ...
# wobei die Variablen immer die Reihenfolge
# a, b, masse, abweichung einhalten.
# Sie sind aufgebaut nach
# <variable1>: (<startwert>, <minimum>, <maximum>)
Z_GAUSS_VARIABLES = (
	(-1, -5, 0),(400, 300, 500),	# Linearer anteil (a*x+b)
	(90, 75, 105), (2, 1, 5)  # Gaußteil Masse, Standardabweichung
)
Z_BREIT_WIGNER_VARIABLES = (
	(-1, -5, 0),(400, 200, 400),
	(90, 75, 105), (4, 2, 10)  # Breit-Wigner Masse, Breite
)
Z_VOIGTIAN_VARIABLES = (
	(-1, -5, 0),(400, 300, 400),
	(90, 75, 105), (2,0.5,5), (1, 0.1, 1.5) # Masse, Breite, Standardabw
)
Z_CRYSTALBALL_VARIABLES = (
	(-1, -5, 0),(400, 300, 400),
	(90, 75, 105), (5,0,10), (100, 0, 200), (100, 0, 200) # masse, sigma, alpha, n
)

J_PSI_GAUSS_VARIABLES = (
	(1, 0, 10),(500, 100, 1000),	# Linearer anteil (a*x+b)
	(3, 2, 4), (0.1, 0.001, 0.5)  # Gaußteil Masse, Standardabweichung
)
J_PSI_BREIT_WIGNER_VARIABLES = (
	(1, 0, 10),(500, 100, 1000),
	(3, 2, 4), (0.1, 0.001, 1)  # Breit-Wigner Masse, Breite
)
J_PSI_VOIGTIAN_VARIABLES = (
	(1, 0, 10),(500, 100, 1000),
	(3, 2, 4), (0.1, 0.01,1), (0.1, 0.001, 0.5) # Masse, Breite, Standardabw
)



def _transformCoordinates(massList, max_m, min_m, nBins=-1, withOverflow = True):
	# in ein Binning transformieren.
	# wenn nBins -1 ist, ist die Anzahl der Bins max_m - min_m
	if nBins == -1: nBins = max_m-min_m
	# Brauchen wir Overflow-Bins?
	# (in Matplotlib sind das erste und das letzte Bin nicht gezeichnete Overflow-Bins)
	step = float(max_m-min_m)/nBins
	if withOverflow:
		x = np.array([(x-1)*step+min_m for x in range(nBins+2)])
		y = np.array([0 for x in range(nBins+2)])
	else:
		x = np.array([(x)*step+min_m for x in range(nBins)])
		y = np.array([0 for x in range(nBins)])
		
	for i in massList:
		Bin = math.floor((i-min_m)/step)
		if withOverflow:
			Bin += 1
		
		try:
			y[Bin] += 1
		except IndexError:
			print "Adde zu Bin %i (habe aber nur Länge %i)"%(Bin, len(y))
			print "Masse %.3f, Minimum %i"%(i, min_m)
			return

	return x,y

def __defaultFitter(massList, roofitfunc, roofitfuncargs,
					plotfitfunc, roo_x, p0,
					outfile, result_blank, min_m, max_m, nBins,
					plt, plotcolor, plotlabel):
	"""
	__defaultFitter - Standard-Fitter

	Funktion, die von allen drei Fitfunktionen (fitGauss,
	fitBreitWigner, fitVoigtian) aufgerufen wird.
	Sie importiert die massList in ein RooFit-Datenformat,
	addiert einen linearen Fit (Hintergrund) zur roofitfunc
	und fittet das anschließend an die Messdaten. Außerdem
	wird das Ergebnis geplottet auf das gegebene Plotobjekt
	und ein ausgewählter Text mit bestimmtem Format (siehe Code)
	wird in die Ausgabedatei des Fits sowie in die Konsole geschrieben.
	Die roofitfuncargs sind hier die Paramter der roofitfunc - wichtig
	zum Schreiben der Ausgabe - es wird hier ein Tupel mit zwei Argumenten
	erwartet.
	"""
	# Binning generieren
	# Under- und Overflow-Bins mitgenerieren - sie existieren im TH1F auch.
	x, y = _transformCoordinates(massList, max_m, min_m, nBins, withOverflow=False)

	# Fehler in den Werten berechnen
	yerrors = np.sqrt(y)

	th = TH1F("invMass", "Invariant Mass", nBins, min_m, max_m)

	# Bininhalte des Histogramms einfügen (Werte und Fehler)
	# Nulltes Bin ist underflow-Bin - einfach nicht beachten.
	for i in range(len(x)):
		th.SetBinContent(i+1, y[i])
		th.SetBinError(i+1, yerrors[i])
		
	# Die Daten in ein RooFit-Datenformat umwandeln
	data = RooDataHist("data", "data", RooArgList(roo_x), th)

	# die Hintergrund-Wahrscheinlichkeitsverteilung
	a = RooRealVar("a", "a", p0[0][0], p0[0][1], p0[0][2])
	b = RooRealVar("b", "b", p0[1][0], p0[1][1], p0[1][2])
	linear = RooPolynomial("background", "Background pdf", roo_x, RooArgList(b, a), 0)

	# Die Anzahlen der Events Hintergrund und Signal
	signalN = RooRealVar("sigN", "Signal number", 10000, 120000)
	backgroundN = RooRealVar("bkgN", "Background number", 500, 80000)

	signalFrac = RooRealVar("sigFrac", "Signal fraction", .8, 0.0, 1.0)

	# Die finale Wahrscheinlichkeitsverteilung und der Fit
	model = RooAddPdf("model", "model", RooArgList(roofitfunc, linear),
					  RooArgList(signalFrac))
					  #RooArgList(signalN, backgroundN))

	print "Fehler im Signal:", signalFrac.getError()
	signalN.setVal(len(massList)*signalFrac.getVal())
	signalN.setError(len(massList)*signalFrac.getError())
	backgroundN.setVal(len(massList)*(1-signalFrac.getVal()))
	backgroundN.setError(len(massList)*signalFrac.getError())

	model.fitTo(data) #, RooFit.Extended())

	# Ergebnistext
	argslist = [backgroundN, a, b, signalN]
	for i in roofitfuncargs: argslist.append(i)
	values = []
	for i in argslist:
		values.append(i.getVal())
		values.append(i.getError())

	result = result_blank%tuple(values)

	outfile.write(result+"\r\n")
	print result

	# Gegebenenfalls plotten
	if not plt == 0:

		# unser x wird 10 mal dichter als die Bins,
		# damit die Kurve schön glatt wird
		# (1. Ableitung einigermaßen stetig)
		nBins *= 10
		step = float(max_m-min_m)/nBins
		plot_x = np.array([x*step+min_m for x in range(nBins)])
		
		plt.plot(plot_x, plotfitfunc(plot_x, argslist), color=plotcolor, label=plotlabel)
	

def fitGauss(massList, outfile, p0, max_m, min_m, nBins=10, plt = 0):
	"""
	fitGauss - Fittet einen linearen Hintergrund auf eine Gaußfunktion

	Führt einen Gaußfit mit einem linearen Hintergrund
	zu den gegebenen Messdaten durch. Plottet den Fit ebenfalls
	auf ein gegebenes matplotlib.pyplot-Objekt.
	Parameter:
	massList: Liste von Elementen, die gemessene Massen enthalten
	outfile: Datei in die der Fitoutput geschrieben wird
	p0: startparameter für die Variablen. Sind als Konstanten definiert
		und können aus diesem Paket importiert werden.
	min_m, max_m: minimale bzw. maximale Masse des Fits
	scalefactor: Körnung des Fits.
	min_m, max_m und der scalefactor ergeben das Binning des Histogramms
	"""

	# Der gültige Wert der x-Variable geht vom Minimum des Diagrams bis zum Maximum
	roo_x = RooRealVar("x", "x", min_m, max_m)
	# Die Annahme für den Mittelwert des Gauss liegt in der Hälfte des Spektrums
	mean = RooRealVar("mean", "mean of gaussian", p0[2][0], p0[2][1], p0[2][2])
	# Standardabweichung
	sigma = RooRealVar("sigma", "width of gaussian", p0[3][0], p0[3][1], p0[3][2])

	gauss = RooGaussian("gauss", "gaussian PDF", roo_x, mean, sigma)

	# Gauß als Python-Funktion zum matplotlib.plotten schreiben
	step = float(max_m-min_m)/nBins
	def gaussFunc(x, args):
		bkg, a,b, sig, mean,sigma = args
		return a.getVal()*x + b.getVal() + \
				sig.getVal()*step/(sigma.getVal()*np.sqrt(2*np.pi))* \
				np.exp(-((x-mean.getVal())**2/(2*sigma.getVal()**2)))

	result = """--- Gaußfit ---
Kurve gefittet. Parameter: a*x + b + Gauss.
#Hintergrund: %i err: %.5f, a: %.5f err: %.5f, b: %.5f err: %.5f
#Prozess: %i err: %.5f, mean: %.5f err: %.5f, sigma: %.5f err: %.5f"""

	__defaultFitter(massList, gauss, (mean, sigma), gaussFunc, roo_x, p0,
					outfile, result, min_m, max_m, nBins,
					plt, "g", "Gaussfit")

def fitBreitWigner(massList, outfile, p0, max_m, min_m, nBins=10, plt = 0):

	step = float(max_m-min_m)/nBins
	def breitWigner(x, params):
		# params = (ax+b + breit wiegner)
		# breit wiegner = m, Max gamma
		bkg, a, b, sig, m, width = params
		m = m.getVal()
		width = width.getVal()
		
		return a.getVal()*x + b.getVal() + \
			   sig.getVal()*step/((x-m)**2+0.25*width*width)

	# Der gültige Wert der x-Variable geht vom Minimum des Diagrams bis zum Maximum
	roo_x = RooRealVar("x", "x", min_m, max_m)
	# Die Annahme für den Mittelwert des Gauss liegt in der Hälfte des Spektrums
	mean = RooRealVar("mean", "mean of breit-wigner", p0[2][0], p0[2][1], p0[2][2])
	# Standardabweichung
	width = RooRealVar("width", "width of breit-wigner", p0[3][0], p0[3][1], p0[3][2])

	bw = RooBreitWigner("breitwigner", "Breit-Wigner PDF", roo_x, mean, width)

	result = """--- Breit-Wigner Fit ---
Kurve gefittet. Parameter: a*x + b + Breit-Wigner.
#Hintergrund: %i err: %.5f, a: %.5f err: %.5f, b: %.5f err: %.5f
#Prozess: %i err: %.5f, mean: %.5f err: %.5f, width: %.5f err: %.5f"""

	__defaultFitter(massList, bw, (mean, width), breitWigner, roo_x, p0,
					outfile, result, min_m, max_m, nBins,
					plt, "r", "Breit-Wigner Fit")

def fitVoigtian(massList, outfile, p0, max_m, min_m, nBins = 10, plt = 0):

	voigtianConst = 1.0/np.sqrt(np.arctan2(0.0, -1.0))

	step = float(max_m-min_m)/nBins
	def voigtianFunc(x, params):
		bgk, a, b, sig, mean, s, w = params

		coef = -0.5/(s.getVal()*s.getVal())
		arg = x-mean.getVal()

		# Voigtian
		c = 1.0/(np.sqrt(2.0)*s.getVal())
		a_ = 0.5*c*w.getVal()
		u = c*arg

		# die roofit implementation
		z = [complex(i, a_) for i in u]
		v = special.wofz(z)

		return a.getVal()*x + b.getVal() + \
			   sig.getVal()*step*c*voigtianConst*v.real

	# RooFit Voigtian + params

	roo_x = RooRealVar("x", "x", min_m, max_m)
	mean = RooRealVar("mean", "mean of voigtian", p0[2][0], p0[2][1], p0[2][2])
	width = RooRealVar("width", "width of breit-wigner-part of vogitian",
					   p0[3][0], p0[3][1], p0[3][2])
	sigma = RooRealVar("sigma", "sigma of gauss-part of vogitian",
					   p0[4][0], p0[4][1], p0[4][2])

	voigt = RooVoigtian("voigtian", "voigtian p.d.f", roo_x, mean, width, sigma)	

	# Ergebnistext
	result = """--- Voigtianfit ---
Kurve gefittet. Parameter: a*x + b + Voigtian.
#Hintergrund: %i err: %.5f, a: %.5f err: %.5f, b: %.5f err: %.5f
#Prozess: %i err: %.5f, mean: %.5f err: %.5f, width: %.5f err: %.5f, sigma: %.5f err: %.5f"""

	__defaultFitter(massList, voigt, (mean, width, sigma), voigtianFunc, roo_x, p0,
					outfile, result, min_m, max_m, nBins,
					plt, "y", "Voigtian Fit")

def fitCrystalball(massList, outfile, p0, max_m, min_m, nBins = 10, plt = 0):

	step = float(max_m-min_m)/nBins
	def crystalballFunc(x, params):
		bgk, a, b, sig, mean, s, alpha, n = params

		n = n.getVal()
		alpha = alpha.getVal()
		t = (x-mean.getVal())/s.getVal()
		if alpha < 0: t = -t

		alpha = math.fabs(alpha)

		l = []
		for i in t:
			if i >= -alpha:
				e = np.exp(-0.5*i*i)
			else:
				a_ = (n/alpha)**n*np.exp(-0.5*alpha*alpha)
				b_ = n/alpha-alpha
				print "n:", n
				e = a_/((b_-t)**n)
			l.append(e)

		return a.getVal()*x + b.getVal() + np.array(e)

	roo_x = RooRealVar("x", "x", min_m, max_m)
	mean = RooRealVar("mean", "mean of crystal ball", p0[2][0], p0[2][1], p0[2][2])
	sigma = RooRealVar("sigma", "sigma of gauss-part of crystal ball",
					   p0[3][0], p0[3][1], p0[3][2])
	alpha = RooRealVar("alpha", "alpha of crystal ball", p0[4][0], p0[4][1], p0[4][2])
	n = RooRealVar("n", "n of crystal ball", p0[5][0], p0[5][1], p0[5][2])

	cb = RooCBShape("cb", "Crystal Ball p.d.f", roo_x, mean, sigma, alpha, n)	

	# Ergebnistext
	result = """--- Crystal Ball fit ---
Kurve gefittet. Parameter: a*x + b + Crystal ball.
#Hintergrund: %i err: %.5f, a: %.5f err: %.5f, b: %.5f err: %.5f
#Prozess: %i err: %.5f, mean: %.5f err: %.5f, sigma: %.5f err: %.5f, alpha: %.5f err: %.5f, n: %.5f err: %.5f"""

	__defaultFitter(massList, cb, (mean, sigma, alpha, n), crystalballFunc, roo_x, p0,
					outfile, result, min_m, max_m, nBins,
					plt, "c", "Crystal-Ball Fit")


	
