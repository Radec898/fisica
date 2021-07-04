#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Módulo controlador del brazo articulado.
# --------------------------------------------------------
# autor: Gabriel Montoiro Varela
# ver:   0.1

import numpy as np
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
            [0, 0, 1, self.d[1], 0, 0] ], dtype=float)

        # matriz del elemento terminal en posición 0
        self.M = np.array([                     \
            [1, 0, 0, 0],                       \
            [0, 1, 0, self.d[1]],               \
            [0, 0, 1, self.d[0] - self.d[2]/2], \
            [0, 0, 0, 1]], dtype=float)

    def get_terminal_pos(self):
        """Devuelve las coordenadas del elemento terminal aplicando cinemática directa."""
        return f2.CinematicaDirectaS(self.M, self.omega.T, self.theta).dot(np.array([0., 0., 0., 1.]))
    
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
        s += "\t - Articulaciones [theta]:\t" + str(self.theta) + "\n"
        s += "\t - Elem. terminal [x,y,z]:\t" + str(self.get_terminal_pos()[:3]) + "\n"
        return s

