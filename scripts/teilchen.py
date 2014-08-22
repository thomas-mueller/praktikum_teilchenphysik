
# -*- coding: utf-8 -*-

import math

"""
teilchen.py

Das Teilchen-Modul.

Autor: Heinrich Blatt
Datum: Dezember 2013

Dieses Modul stellt die Klasse Teilchen zur verfügung. Die Software ist so
ausgelegt, dass man keine weiteren Funktionen bzw. Eigenschaften als diese
Klasse benötigt.
"""

h = 4.135e-15 # eV*s
c = 2.997e8   # m/s

class Teilchen:
	"""
	Diese Klasse beinhaltet für den
	Versuch alle relevanten Daten, die ein Teilchen haben kann. Diese werden mit
	dem Konstruktor übergeben und mit allen anderen Funktionen wird nur noch darauf
	zugegriffen.
	"""
	
	def __init__(self, pt = 0, theta = 0, phi = 0, m = 0, q = 0,
				 numChambers = 0, numPixelhits = 0,
				 numStriphits = 0, chi2DivNDOF = 0,
				 pfIso04=0, eventNum=0, runNum = 0, lumiNum = 0,
				 nVertices = 0, nTracks = 0, entfVertexX = 0,
				 entfVertexY = 0, entfVertexZ = 0):
		"""
		Initialisiert ein neues Teilchen

		pt - Transversaler Impuls des Teilchens
		theta, phi - Raumwinkel des Teilchens (theta = 0 bei Neutrinos)
		m - Masse des Teilchens
		q - Ladung
		numChambers - Anzahl der Kammern, in denen das Teilchen gemessen wurde
		numPixelHits, num StripHits - Anzahl der Pixel / Spurenkammern, in
									  denen ein Hit gemessen wurde
		chi2DivNDOF - Die Spurqualität (s. Arbeit) des Teilchens
		pfIso04 - Isolationsfaktor
		"""

		#pt = float(pt) / (math.cos(math.pi/2+float(theta)))
		
		# 4rer-Impuls als Klassenvariable definieren
		self._theta = float(theta)
		self._phi = float(phi)
		self._pt = float(pt)
		self._m = float(m)

		#self._entfVertex = math.sqrt(float(entfVertexX)*float(entfVertexX)+
		#							 float(entfVertexY)*float(entfVertexY)+
		#							 float(entfVertexZ)*float(entfVertexZ))
		self._entfVertex = float(entfVertexZ)
		 
		if self._theta != 0:
			self._eta = -math.log(math.tan(self._theta/2))
		else:
			self._eta = 0

		# Weitere Eigenschaften als Klassenvariablen definieren
		self._q = int(q)
		self._numChambers = int(numChambers)
		self._numPixelhits = int(numPixelhits)
		self._numStriphits = int(numStriphits)
		self._chi2DivnDOF = float(chi2DivNDOF)

		self._nVertices = int(nVertices)
		self._nTracks = int(nTracks)

		self._pfIso04 = float(pfIso04)
		self._eventNum = eventNum
		self._runNum = runNum
		self._lumiNum = lumiNum

		# Impuls in karthesichen Koordinaten berechnen und als
		# Klassenvariablen speichern
		self._px = self._pt*math.cos(self._phi)
		self._py = self._pt*math.sin(self._phi)

		if self._theta != 0:
			self._pz = self._pt/math.tan(self._theta)
		else:
			self._pz = 0

		# Energie des Teilchens berechnen und speichern
		E = self._m*self._m
		E += self._px*self._px + self._py*self._py + self._pz*self._pz
		self._E = math.sqrt(E)

	def __str__(self):
		"""Eigenschaften des Teilchens als String"""
		return """
Ich bin ein Teilchen.
Mein transversaler Impuls ist %.5f.
Dieser geht in die Raumwinkel Theta = %.5f und Phi = %.5f.
Meine Ladung beträgt %i.
Ich wurde in %i Kammern gemessen, dabei habe ich Hits in %i Pixel- und
in %i Streifendetektoren hinterlassen.
Mein Isolationsfaktor beträgt %.5f.
Meine Energie beträgt %.5f GeV."""%(self._pt, self._theta, self._phi,
									self._q, self._numChambers,
									self._numPixelhits, self._numStriphits,
									self._pfIso04, self._E)

	def getEntfVertex(self):
		"""Entfernung vom Vertex in m"""
		return self._entfVertex*h*c*1e9

	def evtPart(self):
		return self._eventNum + " " + self._runNum + " " + self._lumiNum

	def nVertices(self):
		return self._nVertices

	def nTracks(self):
		return self._nTracks

	def m(self):
		"""Ruhemasse des Teilchens"""
		return self._m
	def E(self):
		"""Energie des Teilchens"""
		return self._E
	def px(self):
		"""x-Impuls des Teilchens"""
		return self._px
	def py(self):
		"""y-Impuls des Teilchens"""
		return self._py
	def pz(self):
		"""z-Impuls des Teilchens"""
		return self._pz
	def p(self):
		"""Gesamtimpuls des Teilchens"""
		return math.sqrt(self._px*self._px+self._py*self._py+self._pz*self._pz)
	def eta(self):
		"""Pseudorapidität des Teilchens"""
		return self._eta

	def pt(self):
		"""Pt of Muon"""
		return self._pt
	def theta(self):
		"""Rapidity of Muon"""
		return self._theta
	def phi(self):
		"""Agle Phi of Muon"""
		return self.phi
	def charge(self):
		"""Ladung des Myons"""
		return self._q
	def numChambers(self):
		"""Anzahl der Kammern die von dem Teilchen getroffen wurden"""
		return self._numChambers
	def numPixelHits(self):
		"""Anzahl an Hits im Pixeldetektor"""
		return self._numPixelhits
	def numStripHits(self):
		"""Anzahl der Hits im Streifendetektor"""
		return self._numStriphits
	def chi2nDOF(self):
		"""Wert von Chi^2/nDOF"""
		return self._chi2DivnDOF

	def isolationFactor(self):
		"""Isolationsfaktor"""
		return self._pfIso04/self._pt
	
	def deltaR(self, other):
		"""DeltaR des Teilchens"""
		dEta = self._eta - other._eta
		dPhi = self._phi - other._phi
		if dPhi < -math.pi: dPhi += 2*math.pi
		if dPhi > math.pi: dPhi -= 2*math.pi
		return math.sqrt(dEta*dEta+dPhi*dPhi)

	def invariantMass(self, other):
		"""Gibt die invariante Masse dieses und des other-Teilchens zurück"""
		m = (other._E+self._E)*(other._E+self._E)
		m -= (other._px+self._px)*(other._px+self._px)
		m -= (other._py+self._py)*(other._py+self._py)
		m -= (other._pz+self._pz)*(other._pz+self._pz)
		return math.sqrt(m)

	def transverseInvariantMass(self, other):
		"""Gibt die invariante Masse dieses und des other-Teilchens zurück"""
		m = (other._E_T()+self._E_T())*(other._E_T()+self._E_T())
		m -= (other._px+self._px)*(other._px+self._px)
		m -= (other._py+self._py)*(other._py+self._py)
		return math.sqrt(m)

	def pT(self):
		return math.sqrt(self._px*self._px+self._py*self._py)

	def _E_T(self):
		p = math.sqrt(self._px*self._px+self._py*self._py+self._pz*self._pz)
		return self.pT()/p*self._E
