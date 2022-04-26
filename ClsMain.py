# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 20:32:44 2017

@author: DARLINGTON MENSAH
         mdarlingtonm@gmail.com
"""
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
import ClsErrorValidation as ev
import time
from scipy.spatial import cKDTree as KDTree


class Main:
    """
    Main class for performing the kriging and surface mass balance
    with the necessary data.
    
    NOTE: For the calculation of z, use the ellipsoidal elevation and not the goidal elevation
    IF not, convert from goidal to ellipsoidal elevation first
    """

    def __init__(self):
        """
        Instantiation of the class and initialisation of class variables
        """
        self._ev = ev

    def openfile(self):
        """
        Function for obtaining a data file path
        Returns a string of file path
        """
        root = tk.Tk()
        root.withdraw()
        filepath = filedialog.askopenfilenames(parent=root, title='Choose a file')
        return filepath

    def importfile(self, filepath):
        """
        Function for importing the prediction grid data file to be used for
        performing kriging interpolation
        Returns:
            Array of prediction data (x,y)
        """
        data = pd.DataFrame()
        for file in filepath:
            df = pd.read_excel(file)
            data = data.append(df)
        data = data.iloc[:, :3]
        data = data[~np.isnan(data).any(axis=1)]
        data = np.array(data, dtype=np.float)

        t = KDTree(data)
        mask = np.ones(data.shape[:1], bool)
        idx = 0
        nxt = 1
        while nxt:
            mask[t.query_ball_point(data[idx], 2)] = False
            nxt = mask[idx:].argmax()
            mask[idx] = True
            idx += nxt
        return data[mask]

    def errorvalidation(self, semi_filepath, krige_filepath, prediction_filepath):
        return self._ev.ErrorValidation().blanking(self.importfile(semi_filepath), 
                                        self.importfile(krige_filepath), self.importfile(prediction_filepath), 0.1, 280)

# pylint: disable=C0103
start_time = time.time()
main = Main()
interpolate = main.errorvalidation(main.openfile(), main.openfile(), main.openfile())

print("---- %s seconds----" % (time.time() - start_time))   
