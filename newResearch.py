import numpy as np
import matplotlib.pyplot as plt
import skgstat as skg
import math


"""
@author: Iván Pérez
         ivan.perezd314@gmail.com
"""


##########################################################################################
###############################       CONSTANTS       ####################################
##########################################################################################

###### GRID SQUARE LENGTH ######
GRID_DIST = 50

###### STARTING GRID VALUES (HANSBREEN) WGS84 ######
#xStartingValue, xFinalValue = 510425.000, 519025.000
#yStartingValue, yFinalValue = 8547775.000, 8564175.000

# Better processing
xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
rows = (yFinalValue-yStartingValue)/GRID_DIST
cols = (xFinalValue-xStartingValue)/GRID_DIST
nTuples = (rows + 1) * (cols + 1)
nTuplesTest = 1000



###### SIGN FUNCTION ######


def sign(x): return math.copysign(1, x)


##### ZERO CONTOUR POINTS #####
zeroContourPoints = []


##########################################################################################
###############################         GRIDS         ####################################
##########################################################################################

# Initial x-axis value: 510425
# Final x-axis value: 519025
# Initial y-axis value: 8547775
# Final y-axis value: 8564175

xDemGrid = np.linspace(xStartingValue, xFinalValue, num=100)
yDemGrid = np.linspace(yStartingValue, yFinalValue, num=200)
vel = np.zeros((173, 329), dtype=float)
xSpeedGrid, ySpeedGrid = [], []
absSpeed = []  # m/year


coords = np.zeros((int(nTuplesTest), 2), dtype=float)
values = np.zeros((int(nTuples), 1), dtype=float)

rr = []
# Data Storage: save speed values and dem grid from files
with open('datafiles/vel_recortada.dat', 'r') as velGrid:
        i = 0
        for line in velGrid.readlines():
            if(float(line.split()[0]) > xFinalValue or float(line.split()[0]) < xStartingValue or float(line.split()[1]) > yFinalValue or float(line.split()[1]) < yStartingValue):
                continue
            x,y = float(line.split()[0]), float(line.split()[1])
            coords[i] = (x, y)
            xRef = int((x-xStartingValue)/GRID_DIST)
            yRef = int((y-yStartingValue)/GRID_DIST)
            # print(xRef, yRef)
            values[i] = float(line.split()[3])
            rr.append(float(line.split()[3]))
            vel[xRef, yRef] = line.split()[3]
            # print(i)
            i += 1
            if( i == nTuplesTest):
                break

print(values)
print(len(values))  

#test = np.fromiter((vel[(c[0] - xStartingValue)/GRID_DIST, (c[1]-yStartingValue)/GRID_DIST] for c in coords), dtype=float)

v = skg.Variogram(coords, rr,maxlag=0.15, n_lags = 30)

v.plot()
v.distance_difference_plot()
plt.show()