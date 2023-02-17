import enum
import numpy as np

xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
GRID_DIST = 50

'''

This method orders speed file in order to decrease execution time.

'''

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

    '''

    This method calculates first aproximation of Bayesian Kriging (monthly), giving speed values and residues.

    '''
    year = 2005
    priorYear = year
    priorMonth = 1
    month = 1

    for i in range(5,78):

        if(i % 12 == 0): #december
            month = 12
        else:
            month = i % 12 

        if( month == 1):
            priorMonth = 12
        else:
            priorMonth = month -1

        dict = {
            1: "Jan",
            2: "Feb",
            3: "Mar",
            4: "Apr",
            5: "May",
            6: "Jun",
            7: "Jul",
            8: "Aug",
            9: "Sep",
            10:"Oct",
            11:"Nov",
            12:"Dec"
        }

        residuesPath = f'speedsAfterKrig/residues/output_residuos_{i}_{year}_2plot.dat'
        firstPriorPath = 'datafiles/vel_interpolada_DEM_ordenada.dat'
        resultPath = f'speedsAfterKrig/speedModules/speedsAfterBayesianKriging_{year}_{dict.get(month)}.dat'
        try:
            priorPath = f'speedsAfterKrig/speedModules/speedsAfterBayesianKriging_{priorYear}_{dict.get(priorMonth)}.dat'
            if(year != priorYear):
                priorYear +=1
        except:
            pass
        with open(residuesPath, 'r') as plottedRes:
            if (month == 5 and year == 2005):
                with open(firstPriorPath, 'r') as prior:
                    with open(resultPath, 'w') as results:
                        resValues = np.loadtxt(plottedRes)
                        priorValues = np.loadtxt(prior)
                    
                        for i, _ in enumerate(resValues): #filas
                            x = resValues[i][0]
                            y = resValues[i][1]
                            speed = resValues[i][2] + priorValues[i][2]/12
                            results.write(f"{x} {y} {speed} \n")
                print("Done")
            else:
                with open(priorPath, 'r') as prior:
                    with open(resultPath, 'w') as results:
                        resValues = np.loadtxt(plottedRes)
                        priorValues = np.loadtxt(prior)

                        for j, _ in enumerate(resValues): #filas
                            x = resValues[j][0]
                            y = resValues[j][1]
                            
                            speed = resValues[j][2] + priorValues[j][2]
                            if(speed < 0 or speed > 100):
                                print(month)
                                print(year)
                                print(resValues[j][2])
                                print(priorValues[j][2])
                            results.write(f"{x} {y} {speed} \n")
        if(month == 12):
            year +=1
            
                    

def main():

    krigingBayesiano()

if __name__ == '__main__':
    main()