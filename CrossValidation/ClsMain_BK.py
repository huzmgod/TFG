# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 20:32:44 2017

@author: DARLINGTON MENSAH
"""
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
import ClsKriging_BK as kg
import CrossValidation.ClsSemivariogram as sv
import time
from scipy.spatial import cKDTree as KDTree


class Main:
    """
    Main class for performing the kriging and surface mass balance
    with the necessary data.
    """

    def __init__(self):
        """
        Instantiation of the class and initialisation of class variables
        """
        self._kg = kg
        self._sv = sv

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
        Function for importing data file. Removes point less than 10cm from
        previous points.
        Returns:
            Array of prediction data (x,y)
        """
        data = pd.DataFrame()
        for file in filepath:
            df = pd.read_excel(file)
            data = data.append(df)
        data = data.iloc[:, :3]
        data = data[~np.isnan(data).any(axis=1)]
        data = np.array(data, dtype=np.float64)
        t = KDTree(data)
        mask = np.ones(data.shape[:1], bool)
        idx = 0
        nxt = 1
        while nxt:
            mask[t.query_ball_point(data[idx], 0.1)] = False
            nxt = mask[idx:].argmax()
            mask[idx] = True
            idx += nxt
        return data[mask]

    def importfile1(self, filepath):
        """
        Function for importing data file. Removes point less than 10cm from
        previous points.
        Returns:
            Array of prediction data (x,y)
        """
        data = pd.DataFrame()
        for file in filepath:
            df = pd.read_excel(file)
            data = data.append(df)
        data = data.iloc[:, :3]
        data = data[~np.isnan(data).any(axis=1)]
        data = np.array(data, dtype=np.float64)
        return data

    def krige(self, semi_filepath, prediction_filepath):
        """
        Function that returns the kriged data for the study data
        INPUT:
            semi_filepath = filepath for the data used in obtaining the semivariogram
            krige_filepath = filepath of experimental data for performing kriging
            prediction_filepath = filepath of grid data of where kriging will be performed.
        OUTPUT:
            File containing the predicted value of each grid coordenates obtained from the kriging operation
        EXTRA INFO:
            If data to be used for obtaining the semivariogram is the same experimental data,
            then semi_filepath = krige_filepath
        """
        krige_filepath = semi_filepath #Para que solo pregunte una vez
        #var_param = self._sv.Semivariogram(semi_filepath).isotropy(nug)
        # input("Press Enter to continue...")
        return self._kg.Kriging().ordinary(krige_filepath, prediction_filepath)

    def ordinary_Krige(self, semi_filepath, prediction_filepath):
       return self.krige(self.importfile(semi_filepath), self.importfile1(prediction_filepath), 0.1)

# pylint: disable=C0103
start_time = time.time()
main = Main()

interpolate = main.ordinary_Krige(main.openfile(), main.openfile())
print("---- %s seconds----" % (time.time() - start_time))
sample = pd.DataFrame(interpolate)
#save_as_name = filedialog.asksaveasfilename() # Meto el nombre del archivo directamente
writer = pd.ExcelWriter('viscosity_interpolated_2.xlsx', engine='xlsxwriter')
sample.to_excel(writer, sheet_name='Kriging')
writer.save()
