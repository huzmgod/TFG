import numpy as np

xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
GRID_DIST = 50

# with open('vel_interpolada_3.dat', 'r') as prior:
#     with open('vel_interpolada_DEM_ordenada_en_y', 'w') as sortedPrior:
#         priorValues = np.loadtxt(prior)
#         currentY = yFinalValue
#         for i, _ in enumerate(priorValues):
#             currentX = xStartingValue
#             for j, value in enumerate(priorValues[i]):
#                 if(value[1] == currentY and value[0]==currentX):
#                     sortedPrior.write(f"{value[0]} {value[1]} {value[5]}")
#                     currentX += GRID_DIST

def sortArray():
    d= []
    # Sort 2D numpy array by 2nd Column
    with open('vel_interpolada_3.dat', 'r') as prior:
        with open('vel_interpolada_DEM_ordenada_en_y.dat', 'w') as sortedPrior:
            for line in prior.readlines():
                line = line.split('  ')
                d.append([float(line[0]),float(line[1]),float(line[5])])
            d = np.array(d)
            print(d)
            sortedArr = d[d[:,1].argsort()]
            print(sortedArr)
            np.savetxt(sortedPrior, sortedArr, fmt='%.4f')
    
def krigingBayesiano():

    with open('output_residuos_2008_mes1_2plot.dat', 'r') as plottedRes:
        with open('vel_interpolada_DEM_ordenada_en_y.dat', 'r') as prior:
            with open('speedsAfterBayesianKriging_2008_m1.dat', 'w') as results:
                resValues = np.loadtxt(plottedRes)
                priorValues = np.loadtxt(prior)
                priorValues = priorValues[priorValues[:,1].argsort()]
                print(priorValues)
                xCoord = xStartingValue
                yCoord = yFinalValue
                finalRes = []
                finalSpeed = []
                checkRes = 0
                checkSpeed = 0
                for j in range(0,173):
                    for i, _ in enumerate(resValues):
                        if(resValues[i][1] < yCoord):
                            continue
                        elif(resValues[i][0] == xCoord and resValues[i][1] == yCoord):
                            finalRes.append(resValues[i][2])
                            checkRes = 1
                        elif(prior[i][0] == xCoord and prior[i][0] == yCoord):
                            finalSpeed.append(prior[i][2]/12)
                            checkSpeed = 0
                        elif(checkSpeed == 1 and checkRes == 1):
                            checkSpeed = 0
                            checkRes = 0
                            xCoord += GRID_DIST
                            yCoord -= GRID_DIST
                            break

                    results.write(f"{resValues[i][0]} {resValues[i][1]} ")

def main():

    krigingBayesiano()

if __name__ == '__main__':
    main()