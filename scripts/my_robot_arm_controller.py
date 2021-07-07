#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Módulo controlador del brazo articulado.
# --------------------------------------------------------
# autor: Gabriel Montoiro Varela
# ver:   0.1

import numpy as np
import sympy as sp
import fisica2 as f2

class MyRoboticArm():
    """Implementa la geometría y lógica de cálculo del brazo articulado."""
    
    def __init__(self, l1, l2, l3):
        """
        Parámetros
        ----------
        l1, l2, l3 : (float)
            Longitudes de los eslabones
        """
	# vector con las dimensiones de los eslabones
        self.d = np.array([l1, l2, l3], dtype=float)
        
        # coordenadas de las articulaciones
        self.theta = np.zeros(6)
        
        # vectores de giro
        self.omega = np.array([                                 \
            [0, 0, 1, 0, 0, 0],                                 \
            [0, 1, 0, -self.d[0], 0, 0],                        \
            [0, 0, 0, 0, 0, 1],                                 \
            [0, 0, 1, self.d[1], 0, 0],                         \
            [1, 0, 0, 0, self.d[0] - self.d[2]/2, -self.d[1]],  \
            [0, 1, 0, self.d[2]/2 - self.d[0], 0, 0] ], dtype=float)

        # matriz del elemento terminal en posición 0
        self.M = np.array([                     \
            [1, 0, 0, 0],                       \
            [0, 1, 0, self.d[1]],               \
            [0, 0, 1, self.d[0] - self.d[2]/2], \
            [0, 0, 0, 1]], dtype=float)

    def get_terminal_pos(self):
        """Devuelve las coordenadas del elemento terminal aplicando cinemática directa."""
        return f2.CinematicaDirectaS(self.M, self.omega.T, self.theta).dot(np.array([0., 0., 0., 1.]))
   
    def get_M(self):
        """Devuelve la matriz de transformación homogénea en la posición 0 del robot."""
        return np.array(self.M)

    def get_S(self):
        """Devuelve los ejes helicoidales en la posición 0 del robot."""
        return np.array(self.omega)
 
    def get_theta(self):
        return np.array(self.theta)

    def set_theta(self, theta):
        """
        Parámetros
        ----------
        theta : []
            coordenadas de las articulaciones
        """
        self.theta = np.array(theta)

    def set_theta_value(self, i, theta):
        """
        Parámetros
        ----------
        i :  int
            id de la articulación
        theta : float
            valor de la articulación
        """
        self.theta[i] = theta

    def reset(self):        
	"""Reinicia el robot a la posición 0."""
        self.theta = np.zeros(6)

    def __str__(self):
        s = "\t[MyRobotArm]\n"
        s += "\t - Eslabones [L1, L2, L3]:\t" + str(self.d) + "\n"
        s += "\t - Articulaciones [theta]:\t" + str(np.round(self.theta, 2)) + "\n"
        s += "\t - Elem. terminal [x,y,z]:\t" + str(np.round(self.get_terminal_pos()[:3], 2)) + "\n"
        return s

    def get_J(self):
        """Devuelve la matriz Jacobiana para la configuración actual del robot."""

        # ejes de rotación
        w = []
        w.append(np.array([0, 0, 1]))
        w.append(np.array([0, 1, 0]))
        w.append(np.array([0, 0, 0]))
        w.append(np.array([0, 0, 1]))
        w.append(np.array([1, 0, 0]))
        w.append(np.array([0, 0, 0])) 

        # vectores de cada eje al siguiente
        q = []
        q.append(np.array([0, 0, 0]))
        q.append(np.array([0, 0, self.d[0]]))
        q.append(np.array([0, self.d[1], 0]))
        q.append(np.array([0, 0, self.d[2]/2]))
        q.append(np.array([0, 0, 0]))
        q.append(np.array([0, 0, 0]))

        # coordendas de las articulaciones
        t = sp.symbols('t0, t1, t2, t3, t4, t5')

        # matrices de rotación a partir de los ejes w (fórmula de Rodrigues)
        R = []
        for i in range(0, 6, 1):
            wmat = f2.VecToso3(w[i])
            R.append(sp.eye(3) + sp.sin(t[i])*wmat + (1 - sp.cos(t[i]))*(wmat*wmat))

        # aplicamos rotaciones a los vectores q,w para llevarlos a la configuración del robot que queremos
        qs = []
        ws = []
        Ri = R[0]
        qs.append(sp.Matrix(q[0]))
        ws.append(sp.Matrix(w[0]))
        for i in range(1, 6, 1):
            ws.append(Ri*sp.Matrix(w[i]))
            qs.append(Ri*sp.Matrix(q[i]) + qs[i-1])
            Ri = Ri * R[i]

        # calculamos las velocidades lineales, los vectores de giro correspondientes y la matriz Jacobiana
        vs = []
        Ji = []
        i = 0
        vs.append(qs[i].cross(ws[i]))
        Ji.append(ws[i].row_insert(3, vs[i]))
        J = Ji[0]
        for i in range(1, 6, 1):
            
            if i==2: vs.append(sp.Matrix([0, 0, 1]))
            else: vs.append(qs[i].cross(ws[i]))
            
            #vs.append(qs[i].cross(ws[i]))
            Ji.append(ws[i].row_insert(3, vs[i]))
            J = J.col_insert(i, Ji[i])

        return J

