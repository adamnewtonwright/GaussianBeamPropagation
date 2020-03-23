# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:51:54 2020

@author: wrighta
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy import oo


# Input Ray parameter, i.e. height and angle
def ray(y,theta):
    '''
    Parameters
    ----------
    y : float or integer or sympy symbol in meters
        The vertical height of a ray.
    theta : float or integer in radians
        The angle of divergence of the ray.

    Returns
    -------
    mat : 2x1 matrix
        [
        [y],
        [teta]
        ]

    '''
    
    mat = np.array([[y],[theta]])
    return mat

# Ray Transfer Matrix for ideal lens with focal length f
def lens(f):
    '''
    Parameters
    ----------
    f : float or integer or sympy symbol in meters
        Thin lens focal length in meters

    Returns
    -------
    mat : 2x2 matrix
    [
    [   1, 0],
    [-1/f, 1]
    ]

    '''

    mat = np.array([[1,0], [-1/f, 1]])
    return mat

# Ray Transfer Matrix for propagation of distance d
def prop(d):
    '''
    Parameters
    ----------
    d : float or integer or sympy symbol
        Distance light is propagating along the z-axis.

    Returns
    -------
    mat: 2x2 matrix
    [
    [1, d],
    [0, 1]
    ]

    '''
    mat = np.array([[1,d], [0,1]])
    return mat


# multiplying the matrices together. mat1 is the last matrix the light interacts with
def mult(mat1,*argv):
    '''
    Parameters
    ----------
    mat1 : 2x2 ABCD matrix
        Last matrix light interacts with.
    *argv : 2x2 ABCD matrices 
        From left to right, the matrices should be entered such that the leftmost matrix interacts
        with light temporally after the rightmost matrix.

    Returns
    -------
    Mat : 2x2 matrix
        The ABCd matrix describing the whole optical system.

    '''

    Mat = mat1
    for arg in argv:
        Mat = np.dot(Mat, arg)
    return Mat

# Adding Gaussian beam parameters
def Zr(wo, lam):
    '''
    Parameters
    ----------
    wo : float, integer, or symbol
        Beam waist radius in meters.
    lam : float, integer, or symbol
        Wavelength of light in meters.

    Returns
    -------
    zr : float, int, symbols
        Rayleigh range for given beam waist and wavelength.

    '''

    zr = np.pi * wo**2 / lam
    return zr

def W0(zr, lam):
    '''
    Parameters
    ----------
    zr : float, integer, symbol
        Rayleigh range in meters
    lam : float, integer, symbol
        Wavelength of light in meters

    Returns
    -------
    w0 : float, integer, symbol
        Beam waist radius in meters

    '''

    w0 = np.sqrt(lam * zr / np.pi)
    return w0

# Remember, there should be an i in front of zr
# but this complicates the calculations, so we usually just let z = 0
# and don't explicitly deal with the i, but still do the math accordingly
#def q0_func(z,zr):
#    qz = z + zr
#    return qz

def q1_func(z, w0, lam, mat):
    '''
    Parameters
    ----------
    z : float, int, symbol
        Position of the beam waist in meters.
    w0 : float, int, symbol
        Radial waist size in meters (of the embedded Gaussian, i.e. W0/M).
    lam : float, int, symbol
        Wavelength of light in meters.
    mat : float, int, symbol
        The ABCD 2x2 matrix describing the optical system.

    Returns
    -------
    z: float, int, symbol
        Position of the beam waist after the optical system
    zr: float, int, symbol
        Rayleigh range of the beam after the optical system
    '''

    A = mat[0][0]
    B = mat[0][1]
    C = mat[1][0]
    D = mat[1][1]
    zr = Zr(w0, lam)
    real = (A*C*(z**2 + zr**2) + z*(A*D + B*C) + B*D) / (C**2*(z**2 + zr**2) + 2*C*D*z + D**2)
    imag = (zr * (A*D - B*C)) / (C**2*(z**2 + zr**2) + 2*C*D*z + D**2)
    z = real
    zr = imag
    return z, zr
    
def q1_inv_func(z, w0, lam, mat):
    '''
    Parameters
    ----------
    z : float, int, symbol
        Position of the beam waist in meters.
    w0 : float, int, symbol
        Radial waist size in meters (of the embedded Gaussian, i.e. W0/M).
    lam : float, int, symbol
        Wavelength of light in meters.
    mat : float, int, symbol
        The ABCD 2x2 matrix describing the optical system.

    Returns
    -------
    R : float, int, symbol
        Radius of curvature of the wavefront in meters.
    w : float, int, symbol
        Radius of the beam in meters.

    '''
    A = mat[0][0]
    B = mat[0][1]
    C = mat[1,0]
    D = mat[1][1]
    zr = Zr(w0, lam)
    real = (A*C*(z**2 + zr**2) + z*(A*D + B*C) + B*D) / (A**2*(z**2 + zr**2) + 2*A*B*z + B**2) 
    imag = -zr * (A*D-B*C) / (A**2 *(z**2 + zr**2) + 2*A*B*z + B**2) 
    R = 1/real
    w = (-lam / imag / np.pi)**.5
    return R, w


def plot(func, var, rang = np.arange(0,3,.01)):
    '''
    Parameters
    ----------
    func : Sympy function of one variable
        Sympy function defining the beam width after the last optical element.
    var : sympy variable
        Variable in func that will be plotted.
    rang : numpy array
        Array of the values along the optical axis to be plotted

    Returns
    -------
    plot : matplotlib graph
        Graph of the beam width of var


    '''
    func = sym.lambdify(var, func)
    plt.figure()
    plt.plot(rang, func(rang), color = 'b')
    plt.plot(rang, -func(rang), color = 'b')
    plt.grid()
    plt.xlabel('Optic Axis (m)')
    plt.ylabel('Beam size (m)')
    plt.show()