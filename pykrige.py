from cmath import nan
import matplotlib.pyplot as plt
import sys
sys.path.append(
    "/home/tracert/miniconda3/pkgs/pykrige-1.6.1-py39h3811e60_1/lib/python3.9/site-packages/pykrige")
from ok import OrdinaryKriging
import kriging_tools as kt
import numpy as np

xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
xDemGrid = np.linspace(xStartingValue, xFinalValue, num=173)
yDemGrid = np.linspace(yStartingValue, yFinalValue, num=329)
GRID_DIST = 50

 
with open('updated_stick_position_2plot.dat', 'r') as stickPositions:
    with open('residuos_mensuales.dat', 'r') as residues:
        with open('zero_contour_points_Dem_50_points_random.dat', 'r') as contourZeros:
            zeros = np.loadtxt(contourZeros)
            resMatrix = np.loadtxt(residues)
            month= 5
            equivalentMonth = month
            year = 2005
            for position in stickPositions.readlines():
                
                print(f'MONTH = {month}')
                print(f'YEAR = {year}')
                # print(f'RES={resMatrix}')
                krigData, aux = [], []
                for i, _ in enumerate(zeros):
                    inputData = [zeros[i][0],zeros[i][1],0.0]
                    krigData.append(inputData)
                row = month - 5 
                positions = position.split(',  ,')
                for i in range(0, 16):
                    stick1 = str(positions[i].split(', '))
                    for char in ['[', ']', ' ', "'", "\n"]:
                        stick1 = stick1.replace(char, "")
                    coords = stick1.split(',')
                    xCoord = float(coords[0])
                    yCoord = float(coords[1])
                    valueRes = resMatrix[row][i]
                    print(f'{valueRes} {str(valueRes) =="nan"}')
                    if (str(valueRes) =="nan" or valueRes == nan \
                        or str(xCoord)=="nan" or str(yCoord)=="nan"): #associated to undefined stick speed
                        continue
                    data = [xCoord, yCoord, valueRes]
                    krigData.append(data)
                # print(krigData)
                krigData = np.array(krigData)
                print(krigData)
                
                OK = OrdinaryKriging(
                    krigData[:, 0],
                    krigData[:, 1],
                    krigData[:, 2],
                    variogram_model="gaussian",
                    verbose=False,
                    enable_plotting=False,
                    exact_values = False
                )
                # except ValueError:
                #     month += 1
                #     if(month % 12 == 0):
                #         year +=1
                #     continue

                ###############################################################################
                # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
                # grid of points, on a masked rectangular grid of points, or with arbitrary points.
                # (See OrdinaryKriging.__doc__ for more information.)

                z, ss = OK.execute("grid", xDemGrid, yDemGrid)

                # ###############################################################################
                # # Writes the kriged grid to an ASCII grid file and plot it.

                kt.write_asc_grid(xDemGrid, yDemGrid, z, filename=f"speedsAfterKrig/output_residuos_{month}_{year}.asc")
                # plt.imshow(z)
                # plt.show()

                with open(f'speedsAfterKrig/output_residuos_{month}_{year}.asc', 'r') as kriggedResInDem:
                    with open(f'speedsAfterKrig/output_residuos_{month}_{year}_clean.dat', 'w') as cleanRes:
                        i=0
                        for line in kriggedResInDem.readlines():
                            if(i>6):
                                cleanRes.write(line)
                            i+=1


                with open(f'speedsAfterKrig/output_residuos_{month}_{year}_clean.dat', 'r') as kriggedResInDem:
                    with open(f'speedsAfterKrig/output_residuos_{month}_{year}_2plot.dat', 'w') as plottedRes:
                        residuosDem=np.loadtxt(kriggedResInDem)
                        currentY = yFinalValue
                        for i, _ in enumerate(residuosDem):
                            currentX = xStartingValue
                            for j, resVal in enumerate(residuosDem[i]):
                                plottedRes.write(f"{currentX} {currentY} {resVal} \n")
                                currentX += GRID_DIST
                            currentY -= GRID_DIST

                month += 1
                if((month % 12) -1== 0):
                    year +=1
                    
                    
                
                
