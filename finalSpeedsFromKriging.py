import enum
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

# def sortArray():
#     d= []
#     # Sort 2D numpy array by 2nd Column
#     with open('vel_interpolada_3.dat', 'r') as prior:
#         with open('vel_interpolada_DEM_ordenada_en_y.dat', 'w') as sortedPrior:
#             for line in prior.readlines():
#                 line = line.split('  ')
#                 d.append([float(line[0]),float(line[1]),float(line[5])])
#             d = np.array(d)
#             print(d)
#             sortedArr = d[d[:,1].argsort()]
#             print(sortedArr)
#             np.savetxt(sortedPrior, sortedArr, fmt='%.4f')

def sortArray():
    d= []
    # Sort 2D numpy array by 2nd Column
    with open('vel_interpolada_3.dat', 'r') as prior:
        with open('vel_interpolada_DEM_ordenada.dat', 'w') as sortedPrior:
            for line in prior.readlines():
                line = line.split('  ')
                d.append([float(line[0]),float(line[1]),float(line[5])])
            d = np.array(d)
            yCoord = yFinalValue
            for i in range(0,329):
                xCoord = xStartingValue
                for j in range(0,173):
                    for row, _ in enumerate(d):
                        if(d[row][0]==xCoord and d[row][1]==yCoord):
                            sortedPrior.write(f'{xCoord} {yCoord} {d[row][2]} \n')
                            print(f'{xCoord} {yCoord} {d[row][2]}')
                            break
                    xCoord += GRID_DIST
                yCoord -= GRID_DIST
            

def krigingBayesiano():
    year = 2005
    for month in range(5,78):
        with open(f'speedsAfterKrig/output_residuos_{month}_{year}_2plot.dat', 'r') as plottedRes:
            with open('vel_interpolada_DEM_ordenada.dat', 'r') as prior:
                with open(f'speedsAfterKrig/speedsAfterBayesianKriging_{month}_{year}.dat', 'w') as results:
                    resValues = np.loadtxt(plottedRes)
                    priorValues = np.loadtxt(prior)
                
                    for i, _ in enumerate(resValues): #filas
                        x = resValues[i][0]
                        y = resValues[i][1]
                        speed = resValues[i][2] + priorValues[i][2]
                        results.write(f"{x} {y} {speed} \n")
        if(month % 12 == 0):
            year +=1
                    

def main():
    # sortArray()

    krigingBayesiano()

if __name__ == '__main__':
    main()