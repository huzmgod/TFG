# -*- coding: utf-8 -*-
"""
Created on Thu Dec 7 11:08:55 2017

@author: Darlington Mensah
         mdarlingtonm@gmail.com
"""
import pandas as pd
import numpy as np
from collections import defaultdict
import time
import ClsSemivariogram as sv
import ClsKriging as kg
from ClsModel import Model as fnc
from pylab import plt
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import warnings


class ErrorValidation:
    """
    ErrorValidation class for performing the interpolation error validation.
    INPUT VALUE:
        semi_filepath: Filepath of file containing data (if elevation, add boundary)
        krige_filepath: Filepath of file containing data (if elevation, DO NOT add boundary) 
        prediction_filepath: Filepath of file containing interpolation grid
    OUTPUT VALUE:
        Returns multiple excel file such as minimum distance and error per each blanking radii file,
        prediction data on grid and useful multiple plot. 
    """

    def __init__(self):
        """
        Instantiation of the class
        """        
        pass

    def blanking(self, semi_data, int_data, prediction, nug, R):
        # =============================================================================
        # Defining and Initializing Variables
        # =============================================================================
        predictor = defaultdict(list)
        mindistance = defaultdict(list)
        color = defaultdict(list)
        frequency = []
        mean = []
        std = []
        lag = []
        classed_resi = []
        classed_dist = []
        col2del = []     
        
        covariance = sv.Semivariogram(semi_data).isotropy(nug)
        
        # =============================================================================
        # R = maximum radius between the two closest data points in the dataset.
        # It is obtained by direct observation of the distribution of the dataset.
        # =============================================================================        
        C = 10
        r0 = R / 100
        sep = np.linspace(R / C, R, C)
        blanking_radii = np.hstack((0, r0, sep))
        
        # =============================================================================
        # Blanking data inside defined radius prior to kriging to obtain interpolation 
        # with its error
        # =============================================================================
        for i in range(len(blanking_radii)):
            krige = []
            min_dist = []
            for j in range(len(int_data)):
                unblanked = semi_data[((semi_data[:, :2] - int_data[j, :2])**2).sum(1) > blanking_radii[i]**2]
                min_dist.append(np.min(cdist(unblanked[:, :2], int_data[j, :2][None])))
                krige.append(np.hstack(kg.Kriging().ordinary(covariance, unblanked, int_data[j, :2])))
                print(str(i) + ' ' + str(j))
             
            predictor[i] = np.hstack(krige)
            mindistance[i] = np.hstack(min_dist)
        del krige
        del min_dist
        
        mindistance = pd.DataFrame(mindistance)
        predictor = pd.DataFrame(predictor)
        color = pd.DataFrame(predictor)
        
        predictor = predictor.T.drop_duplicates().T
        mindistance = mindistance.T.drop_duplicates().T
        color = color.T.drop_duplicates().T

        predictor = predictor.apply(pd.Series.drop_duplicates, axis=1)
        mindistance = mindistance.apply(pd.Series.drop_duplicates, axis=1)
        color = color.apply(pd.Series.drop_duplicates, axis=1)
                
        # =============================================================================
        # Get the interpolator error
        # =============================================================================
        error = (np.array(predictor).transpose() - int_data[:, 2].transpose()).transpose()
 
        # =============================================================================
        # Scatter plot of minimum distance between points versus its interpolator error
        # =============================================================================       
        mindistance = np.array(mindistance)
        color = np.array(color)     

        for i in range(len(predictor.columns)):
            for k in range(len(predictor)):
                color[k, i] = blanking_radii[i]

        upper_x = int(np.ceil(np.max(mindistance)/10.0))*10
        lower_y = int(np.floor(np.min(error)/10.0))*10
        upper_y = int(np.ceil(np.max(error)/10.0))*10   
        
        plt.scatter(mindistance, error, c=color)
        plt.xlim(0, upper_x)
        plt.ylim(lower_y, upper_y)
        plt.savefig('Scatter Plot.pdf', fmt='pdf', dpi=200)
        plt.show()

        for i in range(len(predictor.columns)):
            plt.scatter(mindistance[:, i], error[:, i])
            plt.xlim(0, upper_x)
            plt.ylim(lower_y, upper_y)
            plt.savefig('Scatter-'+str(blanking_radii[i])+'.pdf', fmt='pdf', dpi=200)
            plt.show()
        
        vecresi = np.array(error).ravel()
        vecdist= np.array(mindistance).ravel()
        sep = np.linspace(R / C, R, C)
        lags = (sep[1:] + sep[:-1]) / 2
        lags = np.hstack((0, r0, R / (2 * C), lags, 2*lags[-1]-lags[-2]))
        
        count = -1
        for ilag in lags[:-1]:
            count = count + 1
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                frequency.append(np.sum((vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])))
                classed_resi.append(vecresi[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])])
                classed_dist.append(np.average(vecdist[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])]))
                mean.append(np.average(vecresi[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])]))
                std.append(np.std(vecresi[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])]))
        classed_error = pd.DataFrame(classed_resi).transpose()
        
        lag = np.hstack((0, r0, sep))
        iclassed_error = pd.DataFrame(classed_resi).transpose()
        count = -1
        for i in range(len(classed_error.columns)):
            count = count + 1
            if np.count_nonzero(~np.isnan(np.array(classed_error)[:, i])) < 100:
                col2del.append(count)
                iclassed_error.drop(i, axis=1, inplace=True)
        lag = np.delete(lag, col2del)
        mean = np.delete(mean, col2del)
        std = np.delete(std, col2del)
        classed_dist = np.delete(classed_dist, col2del)
        frequency = np.delete(frequency, col2del)
        
        # =============================================================================
        # Scatter plot of the error gropued in classes defined by the vector "lag"
        # =============================================================================
        iclassed_error = np.array(iclassed_error)
        plt.scatter(np.tile(classed_dist, len(iclassed_error)), iclassed_error.flatten())
        plt.plot(classed_dist, mean, 'o', c='k')
        plt.plot(classed_dist, std, 'o', c='r')
        plt.xlabel('Distance')
        plt.ylabel('Error')
        plt.title('DBF and DEF')
        plt.savefig('DBF and DEF.pdf', fmt='pdf', dpi=200)
        plt.show()
        
        # =============================================================================
        # Least square fit function to obtain coefficient of the indeterminate
        # =============================================================================
        k = 4.6
        a = 3
        b = 1
        num = a-classed_dist*(a-b)/R
        nnum = num/np.sum(num)
            
        denum = 1-np.exp(-frequency*k/d_nrow)
        ndenum = denum/np.sum(denum)

        weight = (nnum * ndenum)
        nweight = weight/np.sum(weight)
        nweight = (1- nweight)/np.sum((1-nweight))
        
        # =============================================================================
        # Least square fit function to obtain coefficient of the indeterminate
        # =============================================================================
        paramsMean = curve_fit(fnc.fitfunction, classed_dist, mean, sigma=nweight, bounds=((0, -np.inf, -np.inf, 0),
                                                            (np.inf, np.inf, np.inf, 0.000001)))
        paramsStD = curve_fit(fnc.fitfunction, classed_dist, std, sigma=nweight, bounds=((-np.inf, -np.inf, -np.inf, 0),
                                                          (np.inf, np.inf, np.inf, nug)))
        
        [m1, m2, m3, m4] = paramsMean[0]
        [s1, s2, s3, s4] = paramsStD[0]
        
        classed_dist = np.hstack((0, classed_dist))
        mean = np.hstack((0, mean))
        std = np.hstack((s4, std))
        frequency = np.hstack((0, frequency))
        x_int = np.linspace(0, np.max(classed_dist), 200)
        mean_int = fnc.fitfunction(x_int, m1, m2, m3, m4)
        std_int = fnc.fitfunction(x_int, s1, s2, s3, s4)
        _, ax = plt.subplots()
        plt.plot(x_int, mean_int)
        plt.plot(x_int, std_int)
        plt.plot(classed_dist, mean, 'o', c='k', label=str(m1) + "x^3 " + str(m2) + "x^2 " + str(m3) + "x")
        plt.plot(classed_dist, std, 'x', c='r', label=str(s1) + "x^3 " + str(s2) + "x^2 " + str(s3) + "x " + str(s4))
        for i, txt in enumerate(frequency):
            ax.annotate(txt, (classed_dist[i], mean[i]))
            ax.annotate(txt, (classed_dist[i], std[i]))
        plt.legend(loc=2, fontsize='xx-small', borderaxespad=0.)
        plt.savefig('Validation Fit.pdf', fmt='pdf', dpi=200)
        plt.show()

        # =============================================================================
        # Estimating the interpolated value for the entire grid of points
        # =============================================================================
        inter = []
        dist_min = []
        for k in range(len(prediction)):
            inter.append(kg.Kriging().ordinary(covariance, semi_data, prediction[k, :2]))
            dist_min.append(np.min(cdist(int_data[:, :2], prediction[k, :2][None])))
            print(str(k) + ' ' + str(len(prediction)))
        inter = np.hstack(inter).T
        dist_min = np.hstack(dist_min)
        krige_mindist = pd.DataFrame(np.column_stack((prediction, inter, dist_min)))
        w_krige_mindist = pd.ExcelWriter(str(datetime.datetime.today())+'_Prediction.xlsx', engine='xlsxwriter')
        krige_mindist.to_excel(w_krige_mindist, sheet_name='Prediction')
        w_krige_mindist.save()
