from cmath import nan
import matplotlib.pyplot as plt
import sys

'''
    Write pykrige installation path using the following example to avoid pykrige problems, where USER is your system username
'''
from pykrige.ok import OrdinaryKriging
import pykrige.kriging_tools as kt
import numpy as np


##################   CONSTANTS   ###########################

xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
xDemGrid = np.linspace(xStartingValue, xFinalValue, num=173)
yDemGrid = np.linspace(yStartingValue, yFinalValue, num=329)
GRID_DIST = 50



##################   ORDINARY KRIGING (RESIDUES)   ###########################


with open('updated_stick_position_2plot.dat', 'r') as stickPositions:
    with open('residuos_mensuales.dat', 'r') as residues:
        with open('datafiles/zero_contour_points_Dem_40_points.dat', 'r') as contourZeros:
            zeros = np.loadtxt(contourZeros)
            resMatrix = np.loadtxt(residues)
            month = 5 #starting data month (May)
            equivalentMonth = month
            year = 2005
            for position in stickPositions.readlines():
                plot = False
                nlags = 20
                if(month == 9 and year == 2008): 
                    plot = True
                    nlags = 56
                print(f'MONTH = {month}')
                print(f'YEAR = {year}')
                # print(f'RES={resMatrix}')
                krigData, aux = [], []
                for i, _ in enumerate(zeros):
                    inputData = [zeros[i][0],zeros[i][1],0.0]
                    krigData.append(inputData)
                row = month - 5 
                positions = position.split(',  ,')
                # print(f'POSITIONS = {positions}')
                for i in range(0, 16):
                    stick1 = str(positions[i].split(', '))
                    for char in ['[', ']', ' ', "'", "\n"]:
                        stick1 = stick1.replace(char, "")
                    # print(f'STICK1 = {stick1}')
                    coords = stick1.split(',')
                    if (coords[0] == ''): coords.pop(0)
                    #format deletion, not much importance
                    # print(f'COORDS = {coords}')
                    xCoord = float(coords[0])
                    yCoord = float(coords[1])
                    valueRes = resMatrix[row][i]
                    if (np.isnan(valueRes) \
                        or np.isnan(xCoord)): #associated to undefined stick speed
                        continue
                    data = [xCoord, yCoord, valueRes]
                    print(f'DATA = {data}')
                    krigData.append(data)
                
                krigData = np.array(krigData)
                  
                OK = OrdinaryKriging(
                    krigData[:, 0],
                    krigData[:, 1],
                    krigData[:, 2],
                    variogram_model="spherical",
                    nlags = nlags,
                    verbose = False,
                    enable_plotting = plot,
                    exact_values = False,
                    pseudo_inv = True
                )
                
                ###############################################################################
                # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
                # grid of points, on a masked rectangular grid of points, or with arbitrary points.
                # (See OrdinaryKriging.__doc__ for more information.)

                z, ss = OK.execute("grid", xDemGrid, yDemGrid)

                # ###############################################################################
                # # Writes the kriged grid to an ASCII grid file and plot it.

                kt.write_asc_grid(xDemGrid, yDemGrid, z, filename=f"speedsAfterKrig/fixedResidues/output_residuos_{year}_{month}.asc")
                # plt.imshow(z)
                # plt.show()

                with open(f'speedsAfterKrig/fixedResidues/output_residuos_{year}_{month}.asc', 'r') as kriggedResInDem:
                    with open(f'speedsAfterKrig/fixedResidues/output_residuos_{year}_{month}_clean.dat', 'w') as cleanRes:
                        i=0
                        for line in kriggedResInDem.readlines():
                            if(i > 6):
                                cleanRes.write(line)
                            i+=1


                with open(f'speedsAfterKrig/fixedResidues/output_residuos_{year}_{month}_clean.dat', 'r') as kriggedResInDem:
                    with open(f'speedsAfterKrig/fixedResidues/output_residuos_{year}_{month}_2plot.dat', 'w') as plottedRes:
                        residuosDem=np.loadtxt(kriggedResInDem)
                        currentY = yFinalValue
                        for i, _ in enumerate(residuosDem):
                            currentX = xStartingValue
                            for j, resVal in enumerate(residuosDem[i]):
                                plottedRes.write(f"{currentX} {currentY} {resVal} \n")
                                currentX += GRID_DIST
                            currentY -= GRID_DIST

                month += 1
                if((month % 12) -1 == 0):
                    year +=1
                    month = 1
                    
                    
                
                
