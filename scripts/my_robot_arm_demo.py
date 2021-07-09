#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Módulo controlador del brazo articulado.
# --------------------------------------------------------
# autor: Gabriel Montoiro Varela
# ver:   0.1

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time
from my_robot_arm_controller import MyRoboticArm
from my_robot_arm_node import MyRoboticArmNode

def clear_screen():
    # códigos de escape ANSI para borrar el terminal (ESC [ 2J) y 
    # posicionarse en la primera línea (ESC [ H)
    print "\033[H\033[2J"

def show_header():
    clear_screen()
    print "\t##########################################"
    print "\t#            F I S I C A   II            #"
    print "\t#                                        #"
    print "\t# - Proyecto: Brazo Articulado 6DOF      #"
    print "\t# - Autor: Gabriel Montoiro Varela       #"
    print "\t# - Curso: 2020/21                       #"
    print "\t##########################################"
    print ""

def show_menu():
    show_header()
    print "\t\tR - Resetear posición"
    print "\t\tE - Mostrar estado"
    print "\t\tP - Posicionar brazo"
    print "\t\tG - Gráfica de posición Elemento Terminal"
    print "\t\tM - Elipsoide de Manipulabilidad/Fuerza"
    print "\t\tS - Salir"
    print ""
    
def get_menu_option():
    opc = None
    while opc not in ('R', 'E', 'P', 'G', 'M', 'S'):
        show_menu()
        opc = raw_input("\tOpción > ").upper()
    return opc

def end_option():
    raw_input("\n\tPulsa [Ret] para salir")

def opc_reset(arm, arm_ros_node):
    """
    Parámetros
    ----------
    arm : MyArmController
        instancia del brazo articulado
    arm_ros_node : MyArmControllerNode
        instancia del nodo ROS
    """
    show_header()

    print "\tReseteando el robot a su posición 0..."

    arm.set_theta(np.zeros(6))
    arm_ros_node.publish_state()

    end_option()

def opc_mostrar_estado(arm):
    """
    Parámetros
    ----------
    arm : MyArmController
        instancia del brazo articulado
    """
    show_header()

    print "\tMostrando estado del robot:\n"
    print arm # imprime el estado completo del robot

    end_option()

def opc_posicionar_brazo(arm, arm_ros_node):
    """
    Parámetros
    ----------
    arm : MyArmController
        instancia del brazo articulado
    arm_ros_node : MyArmControllerNode
        instancia del nodo ROS
    """
    show_header()
    
    print "\tPosicionar elemento terminal del robot:\n"
    print "\tIndica el valor de las articulaciones"
    print "\t - Para las rotaciones el rango máximo debe ser [0, 2*PI]"
    print "\t - Para la traslación el rango máximo debe ser [-L3/2, L3/2]\n"

    theta = arm.get_theta() # inicializamos con los valores actuales
    for i in range(6):
        val = raw_input("\tTheta[" + str(i) + "] (actual=" + str(theta[i]) + "): ")
        if len(val)>0:
            theta[i] = float(val)

    print "\n\tReposicionando brazo..."

    # ajuste por las dimensiones en el archivo URDF
    theta[2] *= (0.13/2.5)
    arm.set_theta(theta)
    arm_ros_node.publish_state()

    print "\n\tNueva posición del elemento terminal:\n"
    print arm

    end_option()

def opc_genera_grafica(arm):
    """
    Genera gráficas 3D por cinemética directa de las posiciones del elemento terminal por
    rotación de las articulaciones 1 y 2 y traslación de la articulación 3
    Las rotaciones de la articulación esférica (4, 5, 6) no se consideran

    - Para las rotaciones el rango máximo debe ser [0, 2PI]
    - Para la traslación el rango máximo debe ser [-L3/2, L3/2]

    Parámetros
    ----------
    arm : MyArmController
        instancia del brazo articulado
    """
    show_header()

    print "\tGenerador de gráficas:\n"
    print "\tIndica las articulaciones a mover, sus valores inicial y final, y el paso"
    print "\t - Para las rotaciones el rango máximo debe ser [0, 2*PI]"
    print "\t - Para la traslación el rango máximo debe ser [-L3/2, L3/2]\n"

    # estructura para los rangos de valores [ini, fin, paso] de las 3 articulaciones
    th_range = np.zeros((3,3))

    # obtenemos los rangos y paso de movimiento
    for i in range(3):
        if raw_input("\tQuieres variar theta" + str(i) + " [S/N]? ").upper() == 'S':
            th_range[i,0] = float(raw_input("\tValor inicial: ")) 
            th_range[i,1] = float(raw_input("\tValor final: ")) 
            th_range[i,2] = float(raw_input("\tPaso: ")) 

    # valores para las articulaciones
    if th_range[0,0] != th_range[0,1]:
        th1_values = np.arange(th_range[0,0], th_range[0,1], th_range[0,2])
    else:
        th1_values = [0.0]
    
    if th_range[1,0] != th_range[1,1]:
        th2_values = np.arange(th_range[1,0], th_range[1,1], th_range[1,2])
    else:
        th2_values = [0.0]
    
    if th_range[2,0] != th_range[2,1]:
        th3_values = np.arange(th_range[2,0], th_range[2,1], th_range[2,2])
    else:
        th3_values = [0.0]

    # obtenemos las coordenadas del elemento terminal para cada combinación de las articulaciones
    pos = []
    old_theta = arm.get_theta()
    arm.reset()
    for t1 in th1_values:
        for t2 in th2_values:
            for t3 in th3_values:
                arm.set_theta([t1, t2, t3, 0., 0., 0.]) # establece valores articulaciones
                pos.append(arm.get_terminal_pos())  # obtiene posición del elemento terminal
    arm.set_theta(old_theta)   
 
    # generamos la gráfica
    arr_coord = np.array(pos) # array de coordenadas generadas
    x = arr_coord[:,0] # lista de coordenada x 
    y = arr_coord[:,1] # lista de coordenada y 
    z = arr_coord[:,2] # lista de coordenada z 

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')
    ax.scatter(x, y, z)
    plt.show()

    end_option()

def opc_genera_elipsoides(arm):
    """
    Genera el elipsoide de manipulabilidad para la configuración actual y velocidades en 
    las tres primeras articulaciones

    Parámetros
    ----------
    arm : MyArmController
        instancia del brazo articulado
    """
    show_header()

    print "\tGenerando esfera de velocidadesi/fuerzas de las articulaciones [θ1/τ1, θ2/τ2, θ3/τ3]...\n"

    # esfera de velocidades de las articulaciones
    r = 1
    u_list = np.linspace(0, 2*np.pi, 25) 
    v_list = np.linspace(0, np.pi, 25) 
    vj1 = np.array([])
    vj2 = np.array([])
    vj3 = np.array([])
    for u in u_list:
        for v in v_list:
            vj1 = np.r_[vj1, r*np.cos(v)*np.cos(u)]
            vj2 = np.r_[vj2, r*np.cos(v)*np.sin(u)]
            vj3 = np.r_[vj3, r*np.sin(v)]
    vj1 = np.r_[vj1, vj1]
    vj2 = np.r_[vj2, vj2]
    vj3 = np.r_[vj3, -vj3]

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel(r'${\dot\theta_1}$ / ${\tau_1}$')
    ax.set_ylabel(r'${\dot\theta_2}$ / ${\tau_2}$')
    ax.set_zlabel(r'${\dot\theta_3}$ / ${\tau_3}$')
    ax.scatter(vj1, vj2, vj3)
    plt.show()

    print "\tGenerando elipsoide de manipulabilidad/fuerza...\n"

    # Jacobiana para la configuración actual
    J = arm.get_analytic_J() # Jacobiana analítica
    #J = arm.get_J() # Jacobiana (ejercicios)
    d = arm.get_D() # dimensiones de los eslabones
    th = arm.get_theta() # posiciones actuales de las articulaciones
    t = sp.symbols('t0, t1, t2, t3, t4, t5')
    L = sp.symbols('L1, L2, L3')
    J_actual = np.array(J.subs({L[0]:d[0], L[1]:d[1], L[2]:d[2], \
                t[0]:th[0],t[1]:th[1],t[2]:th[2],t[3]:th[3],t[4]:th[4],t[5]:th[5]})).astype(float)
    J_inv = np.linalg.pinv(J_actual.T)
   
    # aplicamos la Jacobiana a las velocidades de las articulaciones para obtener las 
    # velocidades del elemento terminal
    vee = np.array([]).reshape([0, 6]) # matriz vacía de 6 columnas para las velocidades del EE
    fee = np.array([]).reshape([0, 6]) # matriz vacía de 6 columnas para las fuerzas en el EE
    for i in range(vj1.shape[0]):
        vjoint = np.array([vj1[i], vj2[i], vj3[i], 0, 0, 0]) # vector de velocidades articulaciones
        vl = np.dot(J_actual, vjoint).reshape([1,6]).astype(float) # nuevo vector de velocidades del EE
        fz = np.dot(J_inv, vjoint).reshape([1,6]).astype(float) # nuevo vector de fuerzas en el EE
        vee = np.r_[vee, vl] # añadimos el nuevo vector de velocidades del EE
        fee = np.r_[fee, fz] # añadimos el nuevo vector de fuerzas en el EE
        
    # mostramos el elipsoide    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel(r'${\dot{x_1}}$ / ${\dot{f_1}}$')
    ax.set_ylabel(r'${\dot{x_2}}$ / ${\dot{f_2}}$')
    ax.set_zlabel(r'${\dot{x_3}}$ / ${\dot{f_3}}$')
    ax.scatter(vee[:,0], vee[:,1], vee[:,2], label='manipulabilidad')
    ax.scatter(fee[:,0], fee[:,1], fee[:,2], label='fuerza')
    plt.legend()
    plt.show()

    end_option()

def main():
    """Programa principal."""

    # creamos una instancia del brazo articulado a partir de las dimensiones de los eslabones
    arm = MyRoboticArm(10, 10, 10)

    # creamos el nodo ROS que publicará el estado del brazo
    arm_ros_node = MyRoboticArmNode(arm)

    time.sleep(1) # pequeña pausa para que se inicie el nodo ROS 
    arm_ros_node.publish_state()
    
    # bucle principal
    opc = None
    while opc != 'S':
        opc = get_menu_option()

        if opc == 'R':
            opc_reset(arm, arm_ros_node)

        if opc == 'E':
            opc_mostrar_estado(arm)

        elif opc == 'P':
            opc_posicionar_brazo(arm, arm_ros_node)     

        elif opc == 'G':
            opc_genera_grafica(arm)     

        elif opc == 'M':
            opc_genera_elipsoides(arm) 

if __name__ == '__main__':
    main()

