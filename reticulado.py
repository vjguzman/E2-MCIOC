# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:00:50 2020

@author: Ale
"""

import numpy as np
from scipy.linalg import solve,inv

class Reticulado(object):
    """Define un reticulado"""
    __NNodosInit__ = 100

    def __init__(self):
        super(Reticulado, self).__init__()
        
        self.xyz = np.zeros((Reticulado.__NNodosInit__,3), dtype=np.double)
        self.Nnodos = 0
        self.barras = []
        self.cargas = {}
        self.restricciones = {}
        self.Ndimensiones = 2
        self.tiene_solucion = False

    def agregar_nodo(self, x, y, z=0):
        if self.Nnodos+1 > Reticulado.__NNodosInit__:
            self.xyz.resize((self.Nnodos+1,3))
        self.xyz[self.Nnodos,:] = [x,y,z]
        self.Nnodos +=1
        if z != 0.:
            self.Ndimensiones = 3
        
    def agregar_barra(self, barra):
        self.barras.append(barra)

    def obtener_coordenada_nodal(self, n): 
        if n >= self.Nnodos:
            return 
        return self.xyz[n, :]

    def calcular_peso_total(self):
        peso = 0.
        for b in self.barras:
            peso += b.calcular_peso(self)
        return peso

    def obtener_nodos(self):
        return self.xyz[0:self.Nnodos,:].copy()

    def obtener_barras(self):
        return self.barras

    def agregar_restriccion(self, nodo, gdl, valor=0.0):
        if nodo not in self.restricciones:
            self.restricciones[nodo] = [[gdl, valor]]
            
        else: 
            self.restricciones[nodo].append([gdl, valor])
            
    def agregar_fuerza(self, nodo, gdl, valor):
        if nodo not in self.cargas:
            self.cargas[nodo] = [[gdl, valor]]
            
        else: 
            self.cargas[nodo].append([gdl, valor])
      

    def ensamblar_sistema(self):
        """Ensambla el sistema de ecuaciones"""
        
        Ngdl = self.Nnodos * self.Ndimensiones

        self.K = np.zeros((Ngdl,Ngdl), dtype=np.double)
        self.f = np.zeros((Ngdl), dtype=np.double)
        self.u = np.zeros((Ngdl), dtype=np.double)
        
        #Implementar
        for b in self.barras:
            ni = b.ni
            nj = b.nj
            d = [2*ni, (2*ni)+1, 2*nj, (2*nj)+1]
            
            for i in range (len(d)):
                p = d[i]
                for j in range(len(d)):
                    q = d[j]

                    ke = b.obtener_rigidez(self)
                    self.K[p,q] += ke[i,j]

                    fe = b.obtener_vector_de_cargas(self)
                    
                self.f[p] += fe[i]
        
        return self.K, self.f 

    def resolver_sistema(self):
        """Resuelve el sistema de ecuaciones.
        La solucion queda guardada en self.u
        """

        # 0 : Aplicar restricciones
        Ngdl = self.Nnodos * self.Ndimensiones
        gdl_libres = np.arange(Ngdl)
        gdl_restringidos = []
        u_otro = np.zeros((Ngdl), dtype=np.double)

        #Identificar gdl_restringidos y llenar u 
        # en valores conocidos.
        #
        # Hint: la funcion numpy.setdiff1d es util
        
        #Pre-llenar el vector u
        for nodo in self.restricciones:
            restriccion = self.restricciones[nodo]
            x = nodo*2
            y = nodo*2 +1
            if len(restriccion) ==2:
                gdl_restringidos.append(x)
                gdl_restringidos.append(y)
                u_otro[x] = restriccion[0][1]
                u_otro[y] = restriccion[1][1]
            else:
                if restriccion[0][0] ==1:
                    gdl_restringidos.append(y)
                    u_otro[y] = restriccion[0][1]
                else:
                    gdl_restringidos.append(x)
                    u_otro[x] = restriccion[0][1]
                
        # con gdl_restringidos encuentro gdl_libres
        gdl_libres = np.setdiff1d(gdl_libres,gdl_restringidos)
        gdl_restringidos = np.array(gdl_restringidos)


        #Agregar cargas nodales a vector de cargas 
        for nodo in self.cargas:
            for carga in self.cargas[nodo]:
                gdl = carga[0]
                valor = carga[1]
                gdl_global = 2*nodo + gdl
                self.f[gdl_global] += valor

        #1 Particionar:
        # K en Kff, Kfc, Kcf y Kcc.
        K = self.K
        Kff = K[np.ix_(gdl_libres,gdl_libres)]
        Kfc = K[np.ix_(gdl_libres,gdl_restringidos)]
        Kcf = Kfc.T

        #f en ff y fc
        ff = []
        F = self.f
        for f in gdl_libres:
            ff.append(F[f])

        #u en uf y uc
        uf = np.setdiff1d(u_otro, gdl_restringidos)
        uc = []
        for valor in u_otro:
            if valor not in uf:
                uc.append(valor)
        uc = np.array(uc)

        # Resolver para obtener uf -->  Kff uf = ff - Kfc*uc
        uf = solve(Kff, ff)
        # porque termino Kfc*uc se va por uc compuesto por ceros
        
        #Asignar uf al vector solucion
        self.u[gdl_libres] = uf

        #Marcar internamente que se tiene solucion
        self.tiene_solucion = True

        return self.u

    def obtener_desplazamiento_nodal(self, n):
        """Entrega desplazamientos en el nodo n como un vector numpy de (2x1) o (3x1)
        """
        dofs = [2*n, 2*n+1]
        return self.u[dofs]

    def recuperar_fuerzas(self):
        """Una vez resuelto el sistema de ecuaciones, se forma un
        vector con todas las fuerzas de las barras. Devuelve un 
        arreglo numpy de (Nbarras x 1)
        """
        
        fuerzas = np.zeros((len(self.barras)), dtype=np.double)
        for i,b in enumerate(self.barras):
            fuerzas[i] = b.obtener_fuerza(self)

        return fuerzas

    def __str__(self):
        s = "nodos:\n"
        for n in range(self.Nnodos):
            s += f"  {n} : ( {self.xyz[n,0]}, {self.xyz[n,1]}, {self.xyz[n,2]}) \n "
        s += "\n\n"

        s += "barras:\n"
        for i, b in enumerate(self.barras):
            n = b.obtener_conectividad()
            s += f" {i} : [ {n[0]} {n[1]} ] \n"
        s += "\n\n"
        
        s += "restricciones:\n"
        for nodo in self.restricciones:
            s += f"{nodo} : {self.restricciones[nodo]}\n"
        s += "\n\n"
        
        s += "cargas:\n"
        for nodo in self.cargas:
            s += f"{nodo} : {self.cargas[nodo]}\n"
        s += "\n\n"

        if self.tiene_solucion:
            s += "desplazamientos:\n"
            if self.Ndimensiones == 2:
                uvw = self.u.reshape((-1,2))
                for n in range(self.Nnodos):
                    s += f"  {n} : ( {uvw[n,0]}, {uvw[n,1]}) \n "
        s += "\n\n"

        if self.tiene_solucion:
            f = self.recuperar_fuerzas()
            s += "fuerzas:\n"
            for b in range(len(self.barras)):
                s += f"  {b} : {f[b]}\n"
        s += "\n"

        return s
