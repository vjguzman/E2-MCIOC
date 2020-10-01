# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:49:08 2020

@author: Ale
"""

import numpy as np

g = 9.81 #kg*m/s^2
phi = 3.14

class Barra(object):

	"""Constructor para una barra"""
	def __init__(self, ni, nj, R, t, E, ρ, σy):
		super(Barra, self).__init__()
		self.ni = ni
		self.nj = nj
		self.R = R
		self.t = t
		self.E = E
		self.ρ = ρ
		self.σy = σy

	def obtener_conectividad(self):
		"""Implementar"""
		return [self.ni,self.nj]

	def calcular_area(self):
		r = self.R - self.t
		area = phi*((self.R**2)-(r**2))
		return area 

	def calcular_largo(self, reticulado):
		"""Devuelve el largo de la barra. 
		xi : Arreglo numpy de dimenson (3,) con coordenadas del nodo i
		xj : Arreglo numpy de dimenson (3,) con coordenadas del nodo j
		"""
		"""saco las coordenadas (x,y,z) de los nodos"""
		n_i = reticulado.obtener_coordenada_nodal(self.ni)
		n_j = reticulado.obtener_coordenada_nodal(self.nj)
		x = abs(n_i[0]-n_j[0])
		y = abs(n_i[1]-n_j[1])
		z = abs(n_i[2]-n_j[2])

		L = (x**2 + y**2)**(0.5)
		return L

	def calcular_peso(self, reticulado):
		"""Devuelve el largo de la barra. 
		xi : Arreglo numpy de dimenson (3,) con coordenadas del nodo i
		xj : Arreglo numpy de dimenson (3,) con coordenadas del nodo j
		"""
		"""Implementar"""
		a = self.calcular_area()
		largo = self.calcular_largo(reticulado)

		peso = (a*largo)*self.ρ*g
		return peso 

	def obtener_rigidez(self, ret):
		"""Devuelve la rigidez ke del elemento. Arreglo numpy de (4x4)
		ret: instancia de objeto tipo reticulado
		"""
		n_i = ret.obtener_coordenada_nodal(self.ni)
		n_j = ret.obtener_coordenada_nodal(self.nj)
		L = self.calcular_largo(ret)

		xi = n_i[0]
		yi = n_i[1]
		xf = n_j[0]
		yf = n_j[1]

		_cos = (xi-xf)/L
		cos = (xf-xi)/L
		_sen = (yi-yf)/L
		sen = (yf-yi)/L

		L = self.calcular_largo(ret)
		A = self.calcular_area()
		k = (self.E * A)/L

		Tθ = np.array([[_cos], [_sen], [cos], [sen]])

		ke = (Tθ @ Tθ.T)*k

		return ke

	def obtener_vector_de_cargas(self, ret):
		"""Devuelve el vector de cargas nodales fe del elemento. Vector numpy de (4x1)
		ret: instancia de objeto tipo reticulado
		"""
		w = self.calcular_peso(ret)

		fe = ((np.array([0.,-1.,0.,-1.]))*(w/2))
		return fe


	def obtener_fuerza(self, ret):
		"""Devuelve la fuerza se que debe resistir la barra. Un escalar tipo double. 
		ret: instancia de objeto tipo reticulado
		"""
		n_i = ret.obtener_coordenada_nodal(self.ni)
		n_j = ret.obtener_coordenada_nodal(self.nj)
		L = self.calcular_largo(ret)
		xi = n_i[0]
		yi = n_i[1]
		xf = n_j[0]
		yf = n_j[1]
		_cos = (xi-xf)/L
		cos = (xf-xi)/L
		_sen = (yi-yf)/L
		sen = (yf-yi)/L

		L = self.calcular_largo(ret)
		A = self.calcular_area()
		k = (self.E * A)/L

		Tθ = np.array([[_cos], [_sen], [cos], [sen]])

		ni = self.ni
		nj = self.nj
		u_e = np.array([ret.u[2*ni], ret.u[2*ni +1], ret.u[2*nj], ret.u[2*nj+1]])
		se = k*(Tθ.T @ u_e)

		return se
	

