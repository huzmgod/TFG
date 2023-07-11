# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 21:00:45 2017

@author: Darlington Mensah
"""
import numpy as np
from math import *
from scipy import integrate

class Model:
    """
    Covariance function model class containing functions for fitting the experimental variogram
    """

    def __init__(self):
        pass
    
    @staticmethod
    
    def posterior_variogram(h):
        mus = 447.969
        mur = 5594.0
        
        sigmas = 82.466
        sigmar = 501.552
        
        c0 = 0.0
        
        # Calculo las 5 integrales para la función de probabilidad conjunta
        
        int1 = lambda s: exp(-(s-mus)**2/(2.0*sigmas**2))
        int2 = lambda r: exp(-(r-mur)**2/(2.0*sigmar**2))
        int3 = lambda s: s*exp(-(s-mus)**2/(2.0*sigmas**2))
        int4 = lambda r: exp(-(r-mur)**2/(2.0*sigmar**2))/r
        int5 = lambda r: exp(-(r-mur)**2/(2.0*sigmar**2))/r**3
        
        integral1 = integrate.quad(int1, 0, 100000)
        integral2 = integrate.quad(int2, 0, 1000000)
        integral3 = integrate.quad(int3, 0, 1000000)
        integral4 = integrate.quad(int4, 0, 1000000)
        integral5 = integrate.quad(int5, 0, 10000)
        
        #print('integral1 = ',integral1)
        #print('integral2 = ',integral2)
        #print('integral3 = ',integral3)
        #print('integral4 = ',integral4)
        #print('integral5 = ',integral5)
        
        # gamma es el variograma posterior, que es función de la distancia h
        
        gamma = np.zeros(len(h))
        for i in range(len(h)):
            if (h[i]<mur):
                gamma[i] = (1/(2*pi*sigmas*sigmar))*(c0*integral1[0]*integral2[0]+integral3[0]*((3*h[i]/2)*integral4[0]-(h[i]**3/2)*integral5[0]))
            else:
                gamma[i] = (1/(2*pi*sigmas*sigmar))*(c0*integral1[0]*integral2[0]+integral3[0]*integral2[0])
        
        return gamma
    
    @staticmethod
    def stable(setlag, setnugget, setsill, setrange, setalpha):
        """
        Stable/Natural covariance function
        Parameters include:
            setlag: a array of the lag of the data
            setnugget: a constant of the nugget value of the data
            setsill: a constant of the sill/variance of the data
            setrange: a constant of the range of the data
            setalpha: a constant which defines the power value (between 1 and 2)
        """   
        return setnugget + setsill*(1-np.exp(-3*(np.array(setlag)/setrange)**setalpha))

    @staticmethod
    def gaussian(setlag, setnugget, setsill, setrange, setalpha=1):
        """
        Gaussian covariance function
        Parameters include:
            setlag: a array of the lag of the data
            setnugget: a constant of the nugget value of the data
            setsill: a constant of the sill/variance of the data
            setrange: a constant of the range of the data
            setalpha: a constant which defines the power value (between 1 and 2)
        """
        del setalpha
        return setnugget + setsill*(1-np.exp(-3*(np.array(setlag)/setrange)**2))

    @staticmethod
    def spherical(setlag, setnugget, setsill, setrange, setalpha=1):
        """
        Spherical covariance function
        Parameters include:
            setlag: a array of the lag of the data
            setnugget: a constant of the nugget value of the data
            setsill: a constant of the sill/variance of the data
            setrange: a constant of the range of the data
            setalpha: a constant which defines the power value (between 1 and 2)
        """
        del setalpha
        setlag = np.array(setlag)
        return np.where(setlag < setrange, setnugget + setsill*(1.5*(setlag/setrange)-
                                     0.5*(setlag/setrange)**3), setnugget + setsill)
        #return np.where(setlag < setrange, setnugget + setsill*(1.5*(np.array(setlag)/setrange)-
        #                             0.5*(np.array(setlag)/setrange)**3), setnugget + setsill)

    @staticmethod
    def exponential(setlag, setnugget, setsill, setrange, setalpha=1):
        """
        Exponential covariance function
        Parameters include:
            setlag: a array of the lag of the data
            setnugget: a constant of the nugget value of the data
            setsill: a constant of the sill/variance of the data
            setrange: a constant of the range of the data
            setalpha: a constant which defines the power value (between 1 and 2)
        """
        del setalpha
        return setnugget + setsill*(1-np.exp(-3*(np.array(setlag)/setrange)))

    @staticmethod
    def fitfunction(x, a, b, c, d):
        """
        3º Polynomial function for fitting mean and StD funcion via least square fit
        """
        return a * x**3 + b * x**2 + c * x + d
