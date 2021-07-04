# -*- coding: utf-8 -*-

# Módulo de funciones para Física II.
# --------------------------------------------------------
# autor: Gabriel Montoiro Varela
# ver:   0.1

import numpy as np

def Rotx(tetha):
    """Crea una matriz de rotación sobre el eje X."""
    c = np.cos(tetha)
    s = np.sin(tetha)
    return np.array(
        [
            [1, 0,  0],
            [0, c, -s],
            [0, s,  c]
        ])

def Roty(tetha):
    """Crea una matriz de rotación sobre el eje Y."""
    c = np.cos(tetha)
    s = np.sin(tetha)
    return np.array(
        [
            [ c, 0, s],
            [ 0, 1, 0],
            [-s, 0, c]
        ])

def Rotz(tetha):
    """Crea una matriz de rotación sobre el eje Z."""
    c = np.cos(tetha)
    s = np.sin(tetha)
    return np.array(
        [
            [c, -s, 0],
            [s,  c, 0],
            [0,  0, 1]
        ])

def norm(v):
    """Normaliza un vector"""
    return v/np.linalg.norm(v)

def Rotv(w, tetha):
    """Crea una matriz de rotación para el eje y ángulo suministrado."""
    w1, w2, w3 = norm(w)
    c = np.cos(tetha)
    s = np.sin(tetha)
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

def MatrixExp3(w, tetha):
    """
    Devuelve la matriz exponencial a partir de un eje w y un ángulo tetha.
    """
    so3 = VecToso3(w)
    return np.identity(3) + np.sin(tetha)*so3 + (1 - np.cos(tetha))*so3.dot(so3)

def MatrixLog3(R):
    if np.array_equal(R, np.identity(3)):
        return .0, None

    if np.trace(R) == -1:
        w1 = (1 / np.sqrt(2 * (1 + R[2,2]))) * np.array([R[0,2], R[1,2], 1 + R[2,2]])
        w2 = (1 / np.sqrt(2 * (1 + R[1,1]))) * np.array([R[0,1], 1 + R[1,1], R[2,1]])
        w3 = (1 / np.sqrt(2 * (1 + R[0,0]))) * np.array([1 + R[0,0], R[1,0], R[2,0]])
        return np.pi, w1, w2, w3

    tetha = np.arccos(0.5 * (np.trace(R) - 1))
    anti_w = (1 / (2 * np.sin(tetha))) * (R - R.T)
    w = so3ToVec(anti_w)
    return tetha, w

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
    v_tetha = se3mat[:3, 3] # v*tetha
    wmat_tetha = se3mat[:3, :3] # [w]*tetha
    w_tetha = so3ToVec(wmat_tetha)
    tetha = np.linalg.norm(w_tetha)

    if tetha < 1E-6:
        # no hay giro
        return np.r_[np.c_[np.identity(3), v_tetha], [[0, 0, 0, 1]]]
    else:
        wmat = wmat_tetha/tetha
        R_wtheta = np.identity(3) + np.sin(tetha)*wmat + (1 - np.cos(tetha))*np.dot(wmat, wmat)
        G_wtetha = np.identity(3)*tetha + (1 - np.cos(tetha))*wmat + \
                    (tetha - np.sin(tetha))*np.dot(wmat, wmat)
        return np.r_[np.c_[R_wtheta, G_wtetha.dot(v_tetha/tetha)], [[0, 0, 0, 1]]]

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

