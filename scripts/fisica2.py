# -*- coding: utf-8 -*-

# Módulo de funciones para Física II.
# --------------------------------------------------------
# autor: Gabriel Montoiro Varela
# ver:   0.1

import numpy as np

def Rotx(theta):
    """Crea una matriz de rotación sobre el eje X."""
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [
            [1, 0,  0],
            [0, c, -s],
            [0, s,  c]
        ])

def Roty(theta):
    """Crea una matriz de rotación sobre el eje Y."""
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [
            [ c, 0, s],
            [ 0, 1, 0],
            [-s, 0, c]
        ])

def Rotz(theta):
    """Crea una matriz de rotación sobre el eje Z."""
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [
            [c, -s, 0],
            [s,  c, 0],
            [0,  0, 1]
        ])

def norm(v):
    """Normaliza un vector"""
    return v/np.linalg.norm(v)

def Rotv(w, theta):
    """Crea una matriz de rotación para el eje y ángulo suministrado."""
    w1, w2, w3 = norm(w)
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [
            [c + w1*w1*(1 - c),     w1*w2*(1 - c) - w3*s,   w1*w3*(1 - c) + w2*s],
            [w1*w2*(1 - c) + w3*s,  c + w2*w2*(1 - c),      w2*w3*(1 - c) - w1*s],
            [w1*w3*(1 - c) - w2*s,  w2*w3*(1 - c) + w1*s,   c + w3*w3*(1 - c)]
        ])

def VecToso3(v):
    """Devuelve la matriz antisimétrica a partir del vector v."""
    return np.array(
        [
            [   0.0, -v[2],  v[1] ],
            [  v[2],   0.0, -v[0] ],
            [ -v[1],  v[0],   0.0 ]
        ])

def so3ToVec(v):
    """Devuelve el vector a partir de la matriz antisimétrica."""
    return np.array([v[2,1], v[0,2], v[1,0]])

def MatrixExp3(w, theta):
    """
    Formula de Rodrigues: Calcula la matriz exponencial de [w]*th a partir de un eje w y un ángulo theta.
    """
    so3 = VecToso3(w)
    return np.identity(3) + np.sin(theta)*so3 + (1 - np.cos(theta))*so3.dot(so3)

def MatrixLog3_parts(R):
    """
    Devuelve el eje unitario y el ángulo de rotación de una matriz SO(3)
    """
    if np.array_equal(R, np.identity(3)):
        return .0, None

    if np.trace(R) == -1:
        w1 = (1 / np.sqrt(2 * (1 + R[2,2]))) * np.array([R[0,2], R[1,2], 1 + R[2,2]])
        w2 = (1 / np.sqrt(2 * (1 + R[1,1]))) * np.array([R[0,1], 1 + R[1,1], R[2,1]])
        w3 = (1 / np.sqrt(2 * (1 + R[0,0]))) * np.array([1 + R[0,0], R[1,0], R[2,0]])
        return np.pi, np,array([w1, w2, w3])

    theta = np.arccos(0.5 * (np.trace(R) - 1))
    anti_w = (1 / (2 * np.sin(theta))) * (R - R.T)
    w = so3ToVec(anti_w)
    return theta, w

def MatrixLog3(R):
    """
    Devuelve la matriz logaritmo de una matriz de rotación
    """
    acosinput = (np.trace(R) - 1) *0.5 
    if np.trace(R) >= 3: return np.zeros((3, 3))
    elif np.trace(R) <= -1:
        if abs(1 + R[2][2])>1.e-6:   omg = (1.0 / np.sqrt(2 * (1 + R[2][2]))) * np.array([R[0][2], R[1][2], 1 + R[2][2]])
        elif abs(1 + R[1][1])>1.e-6: omg = (1.0 / np.sqrt(2 * (1 + R[1][1]))) * np.array([R[0][1], 1 + R[1][1], R[2][1]])
        else: omg = (1.0 / np.sqrt(2 * (1 + R[0][0]))) * np.array([1 + R[0][0], R[1][0], R[2][0]])
        return VecToso3(np.pi * omg)
    else:
        theta = np.arccos(acosinput)
        return (theta*0.5)/np.sin(theta) * (R-np.array(R).T)

def RpToTrans(R, p):
    """
    Convierte una matriz de rotación R y un vector de posición p en una
    matriz homogénea.
    """
    n = R.shape[0]
    H = np.zeros((n + 1, n + 1))
    H[:n, :n] = R
    H[:n, n] = p
    H[n, n] = 1
    return H

def TransToRp(T):
    """
    Convierte una matriz homogénea en una matriz de rotación y un vector de
    posición.
    """
    n = T.shape[0]
    R = T[:n-1, :n-1]
    p = T[:n-1, n-1]
    return R, p

def TransInv(T):
    """
    Inversa de la matriz de transformación homogénea T
    """
    R,p = TransToRp(T)
    Tinv = RpToTrans(R.T, -R.T.dot(p))
    return Tinv

def VecTose3(V):
    """
    Convierte un vector giro v ∈ R6 an una matriz ∈ se(3).
    """
    w = V[:3]
    v = V[3:]
    w_anti = VecToso3(w)
    M = np.zeros((4, 4))
    M[:3, :3] = w_anti
    M[0:3, 3] = v
    return M

def se3ToVec(se3mat):
    """
    Convierte la matriz se(3) en un vector giro ∈ R6
    """
    w = so3ToVec(se3mat[:3, :3])
    v = se3mat[:3, 3]
    return np.concatenate([w, v])

def Adjunta(T):
    """
    Devuelve la adjunta de un matriz de transformación homogénea.
    """
    R = T[:3, :3]
    p = T[:3, 3]
    A = np.zeros((6, 6))
    A[:3, :3] = R
    A[3:, 3:] = R
    A[3:, :3] = VecToso3(p).dot(R)
    return A

def MatrixExp6(se3mat):
    """
    Matriz exponencial de un vector de giro en representación matricial 4x4.
    """
    v_theta = se3mat[:3, 3] # v*theta
    wmat_theta = se3mat[:3, :3] # [w]*theta
    w_theta = so3ToVec(wmat_theta)
    theta = np.linalg.norm(w_theta)

    if theta < 1E-6:
        # no hay giro
        return np.r_[np.c_[np.identity(3), v_theta], [[0, 0, 0, 1]]]
    else:
        wmat = wmat_theta/theta
        R_wtheta = np.identity(3) + np.sin(theta)*wmat + (1 - np.cos(theta))*np.dot(wmat, wmat)
        G_wtheta = np.identity(3)*theta + (1 - np.cos(theta))*wmat + \
                    (theta - np.sin(theta))*np.dot(wmat, wmat)
        return np.r_[np.c_[R_wtheta, G_wtheta.dot(v_theta/theta)], [[0, 0, 0, 1]]]

def MatrixLog6(T):
    """
    Calcula la matriz logaritmo de una matriz de transformación homogénea
    """
    R = T[:3, :3]   # rotación
    p = T[:3, 3]    # traslación

    omgmat = MatrixLog3(R) # coordenadas exponenciales de la matriz de rotación
                           # o sea, un vector de rotación como matriz antisimétrica so3 (3x3)
    if np.array_equal(omgmat, np.zeros((3, 3))): # Si no hay rotación, es una matriz de ceros 
        return np.r_[np.c_[np.zeros((3, 3)),p],[[0, 0, 0, 0]]]
    else:
        omgvec= so3ToVec(omgmat) # expresa la rotación como un vector en la dirección del eje por el ángulo
        omgmat=omgmat/np.linalg.norm(omgvec) # el vector en el eje de rotación normalizado y en forma matricial
        theta = np.linalg.norm(omgvec) # también se puede calcular como np.arccos((np.trace(R)-1)/2.0)
        # a continuación aplicamos la definición que vimos en clase (ver diapositivas)
        invG_theta=np.eye(3)/theta-omgmat*0.5+(1.0/theta-0.5/np.tan(theta*0.5))*np.dot(omgmat,omgmat)
        v=np.dot(invG_theta,p)
        return np.r_[np.c_[omgmat,v],[[0, 0, 0, 0]]]*theta # primero concatena columnas y luego filas    

def CinematicaDirectaS(M, Smatrix, theta):
    """Devuelve la matriz de transformación homogénea del elemento terminal.

    Parámetros
    ----------
    M : np.array([4,4])
        Matriz homogénea del elemento terminal en posición 0 del robot
    Smatrix : np.array([6,n_dof])
        Vectores de giro (columnas)
    theta : np.array(n_dof)
        Coordenadas de las articulaciones
    """
    n_dof = Smatrix.shape[1]

    # Matrices de giro [Si]
    S = np.zeros((n_dof, 4, 4))
    for i in range(n_dof):
        S[i] = VecTose3(Smatrix[:,i])

    # Matrices exponenciales
    E = np.zeros((n_dof, 4, 4))
    for i in range(n_dof):
        E[i] = MatrixExp6(S[i]*theta[i])

    # Transformacion
    T = E[0]
    for i in range(1, n_dof):
        T = T.dot(E[i])

    return T.dot(M)

def TransMatrix(theta, p):
    """Devuelve la matriz de transformación homogénea dada la orientación y la posición.

    Parámetros
    ----------
    theta : np.array(3)
        orientación en los tres ejes (ángulos de Euler)
    p : np.array(3)
        posición (x, y, z)
    """

    # ejes
    eje_x = np.array([1, 0, 0])
    eje_y = np.array([0, 1, 0])
    eje_z = np.array([0, 0, 1])

    # matrices exponenciales por cada eje
    Rx = MatrixExp3(eje_x, theta[0])
    Ry = MatrixExp3(eje_y, theta[1])
    Rz = MatrixExp3(eje_z, theta[2])

    # matriz de rotación
    R = Rz.dot(Ry.dot(Rx))
    print R, p

    # devolvemos la matriz homogénea
    return np.r_[np.c_[R, p], np.array([[0, 0, 0, 1]])]

