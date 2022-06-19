"""
@author: Iván Pérez
         ivan.perezd314@gmail.com
"""


import numpy as np
import math
from pyrsistent import v
from sqlalchemy import null

##########################################################################################
###############################       CONSTANTS       ####################################
##########################################################################################

###### GRID SQUARE LENGTH ######
GRID_DIST = 50

###### STARTING GRID VALUES (HANSBREEN) WGS84 ######
xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
rows = (yFinalValue-yStartingValue)/GRID_DIST
cols = (xFinalValue-xStartingValue)/GRID_DIST

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

xDemGrid = np.linspace(xStartingValue, xFinalValue, num=173)
yDemGrid = np.linspace(yStartingValue, yFinalValue, num=329)
zDemGrid = np.zeros((329, 173), dtype=float)
xSpeedGrid, ySpeedGrid = [], []
absSpeed = []  # m/year

# Data Storage: save speed values and dem grid from files
with open('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat', 'r') as demGrid:
    with open('vel_recortada.dat', 'r') as velGrid:
        with open('cuadricula_superficie_Hansbreen.dat', 'w') as sortedGrid:

            for line in demGrid.readlines():
                value = line.split()
                yIndex = yStartingValue
                y = int((float(value[0]) - xStartingValue)/GRID_DIST)
                # descendent order from greater y to lower y
                x = 328 - int((float(value[1]) - yStartingValue)/GRID_DIST)
                zDemGrid[x, y] = value[2]
            for row, value in enumerate(zDemGrid):

                sortedGrid.write(f"{' '.join(map(str,value))}\n")

            for line in velGrid.readlines():
                if(float(line.split()[0]) > xFinalValue or float(line.split()[0]) < xStartingValue or float(line.split()[1]) > yFinalValue or float(line.split()[1]) < yStartingValue):
                    continue
                xSpeedGrid.append(line.split()[0])
                ySpeedGrid.append(line.split()[1])
                absSpeed.append(line.split()[3])

#### FUNCTION TO INTERPOLATE SPEED GRID INTO DEM GRID ####


def interpolate():
    # Weighted average distance interpolation
    with open('vel_interpolada_3.dat', 'w') as outputSpeed:
        for i, x in enumerate(xDemGrid):
            for _, y in enumerate(yDemGrid):
                coordXdem = float(x)
                coordYdem = float(y)
                distmin = [0.0, 0.0, 0.0]
                speed = [0.0, 0.0, 0.0]
                # vector of closest points
                p = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
                denom = 0.0
                estimatedSpeed = 0.0
                closePoint = False

                for k, _ in enumerate(xSpeedGrid):
                    coordXvel = float(xSpeedGrid[k])
                    coordYvel = float(ySpeedGrid[k])
                    modSpeed = float(absSpeed[k])
                    dist = math.sqrt((coordXdem-coordXvel) **
                                     2+(coordYdem-coordYvel)**2)
                    if(distmin[0] == 0):  # first entry
                        distmin[0] = dist
                        p[0] = (coordXvel, coordYvel)
                        speed[0] = modSpeed
                    elif(dist < distmin[0]):  # getting minimum distance
                        distmin[2] = distmin[1]
                        p[2] = p[1]
                        speed[2] = speed[1]
                        distmin[1] = distmin[0]
                        p[1] = p[0]
                        speed[1] = speed[0]
                        distmin[0] = dist
                        p[0] = (coordXvel, coordYvel)
                        speed[0] = modSpeed
                    elif(distmin[1] == 0):  # second entry
                        distmin[1] = dist
                        p[1] = (coordXvel, coordYvel)
                        speed[1] = modSpeed
                    elif(dist < distmin[1]):  # getting minimum distance
                        distmin[2] = distmin[1]
                        p[2] = p[1]
                        speed[2] = speed[1]
                        distmin[1] = dist
                        p[1] = (coordXvel, coordYvel)
                        speed[1] = modSpeed
                    elif(distmin[2] == 0 or dist < distmin[2]):  # third entry
                        distmin[2] = dist
                        p[2] = (coordXvel, coordYvel)
                        speed[2] = modSpeed

                # Loop that calculates sum of square distances
                for i, dist in enumerate(distmin):
                    # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                    if(dist < 0.5):
                        estimatedSpeed = speed[i]
                        closePoint = True
                        continue
                    else:
                        denom += 1/dist**2

                # if point from speed grid is not so close to point on DEM:
                # there are no close points
                if(closePoint == False):
                    for i, dist in enumerate(distmin):
                        # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                        estimatedSpeed += 1/dist**2/denom*speed[i]

                outputSpeed.write(
                    f"{coordXdem}  {coordYdem}  {p}  {distmin}  {speed}  {estimatedSpeed} \n")


def gradient():
    # Check gradient output
    gradx, grady = np.gradient(zDemGrid)
    with open("gradient_x_dem.dat", 'w') as gradXFile:
        with open("gradient_y_dem.dat", 'w') as gradYFile:
            np.savetxt(gradXFile, gradx, fmt='%.4f')
            np.savetxt(gradYFile, grady, fmt='%.4f')


def speedComponentsDem(gradX, gradY):  # gradX and gradY are .dat files in matrix form
    # This method calculates speed components for every single coordinate in DEM, based on gradient values
    with open("vel_interpolada_3.dat", 'r') as speedValues:
        with open(f"{gradX}", 'r') as gradXFile:
            with open(f"{gradY}", 'r') as gradYFile:
                with open("speed_components.dat", 'w') as speedFile:
                    gradx = np.loadtxt(gradXFile)
                    grady = np.loadtxt(gradYFile)

                    for line in speedValues.readlines():
                        xCoord = float(line.split('  ')[0])
                        yCoord = float(line.split('  ')[1])
                        equalCoords = False
                        for i, _ in enumerate(gradx):
                            for j, _ in enumerate(gradx[i]):
                                xMatrix = xStartingValue + j*GRID_DIST
                                yMatrix = yFinalValue - i*GRID_DIST
                                if (xMatrix == xCoord and yMatrix == yCoord):
                                    equalCoords = True
                                    xValue = gradx[i, j]
                                    yValue = grady[i, j]
                                    xSpeed = 0.0
                                    ySpeed = 0.0
                                    speed = float(line.split('  ')[5])

                                    hypotenuse = math.sqrt(
                                        xValue**2 + yValue**2)

                                    if(xValue == 0 and yValue != 0):
                                        ySpeed = speed*sign(yValue)

                                    elif(yValue == 0 and xValue != 0):
                                        xSpeed = speed*sign(xValue)

                                    elif(xValue != 0 and yValue != 0):
                                        xSpeed = speed*xValue/hypotenuse
                                        ySpeed = speed*yValue/hypotenuse

                                    speedFile.write(
                                        f"{xCoord}  {yCoord}  {xMatrix}  {yMatrix}  {xSpeed}  {ySpeed}\n")
                                    break
                            if(equalCoords == True):
                                break


# calculates gradients for all sticks in first instance
def getGradientInStickCoordinates(filenameX, filenameY):

    ###### STICKS COORDINATES ######

    with open('posiciones_tyczki.dat', 'r') as posicionesEstacas:
        with open(f'{filenameX}', 'r') as gradXFile:
            with open(f'{filenameY}', 'r') as gradYFile:
                with open('gradient_in_stick.dat', 'w') as gradStickFile:
                    gradx = np.loadtxt(gradXFile)
                    grady = np.loadtxt(gradYFile)

                    for line in posicionesEstacas.readlines():
                        distmin = [0.0, 0.0, 0.0]
                        xCoord = float(line.split()[0])
                        yCoord = float(line.split()[1])
                        # vector of closest points
                        p = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
                        denom = 0.0
                        gradValues = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
                        closePoint = False
                        estimatedGradientX, estimatedGradientY = 0.0, 0.0
                        for i, _ in enumerate(gradx):
                            for j, _ in enumerate(gradx[i]):
                                xMatrix = xStartingValue + j*GRID_DIST
                                yMatrix = yFinalValue - i*GRID_DIST
                                dist = math.sqrt(
                                    (xMatrix-xCoord)**2 + (yMatrix-yCoord)**2)

                                if(distmin[0] == 0):  # first entry
                                    distmin[0] = dist
                                    p[0] = (xMatrix, yMatrix)
                                    gradValues[0] = (gradx[i, j], grady[i, j])
                                elif(dist < distmin[0]):  # getting minimum distance
                                    distmin[2] = distmin[1]
                                    p[2] = p[1]
                                    gradValues[2] = gradValues[1]
                                    distmin[1] = distmin[0]
                                    p[1] = p[0]
                                    gradValues[1] = gradValues[0]
                                    distmin[0] = dist
                                    p[0] = (xMatrix, yMatrix)
                                    gradValues[0] = gradx[i, j], grady[i, j]
                                elif(distmin[1] == 0):  # second entry
                                    distmin[1] = dist
                                    p[1] = (xMatrix, yMatrix)
                                    gradValues[1] = gradx[i, j], grady[i, j]
                                elif(dist < distmin[1]):  # getting minimum distance
                                    distmin[2] = distmin[1]
                                    p[2] = p[1]
                                    gradValues[2] = gradValues[1]
                                    distmin[1] = dist
                                    p[1] = (xMatrix, yMatrix)
                                    gradValues[1] = gradx[i, j], grady[i, j]
                                elif(distmin[2] == 0 or dist < distmin[2]):  # third entry
                                    distmin[2] = dist
                                    p[2] = (xMatrix, yMatrix)
                                    gradValues[2] = gradx[i, j], grady[i, j]

                        # Loop that calculates sum of square distances
                        for i, dist in enumerate(distmin):
                            # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                            if(dist < 0.5):
                                estimatedGradientX = gradValues[i][0]
                                estimatedGradientY = gradValues[i][1]
                                closePoint = True
                                continue
                            else:
                                denom += 1/dist**2

                        # if stick point is not so close to point on DEM:
                        # there are no close points
                        if(closePoint == False):
                            for i, dist in enumerate(distmin):
                                # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                                estimatedGradientX += 1 / \
                                    dist**2/denom*gradValues[i][0]
                                estimatedGradientY += 1 / \
                                    dist**2/denom*gradValues[i][1]

                        gradStickFile.write(
                            f"{xCoord}  {yCoord}  {p}  {distmin}  {estimatedGradientX}  {estimatedGradientY} \n")


def getGradientInSpecificStick(gradX, gradY, stickPosition):

    ###### This method calculates weighted average gradient around stick position ######

    with open(f'{gradX}', 'r') as gradXFile:
        with open(f'{gradY}', 'r') as gradYFile:

            gradx = np.loadtxt(gradXFile)
            grady = np.loadtxt(gradYFile)

            distmin = [0.0, 0.0, 0.0]

            xCoord = stickPosition[0]
            yCoord = stickPosition[1]

            # vector of closest points
            p = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
            denom = 0.0
            gradValues = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
            closePoint = False
            estimatedGradientX, estimatedGradientY = 0.0, 0.0

            for i, _ in enumerate(gradx):
                for j, _ in enumerate(gradx[i]):

                    xMatrix = xStartingValue + j*GRID_DIST
                    yMatrix = yFinalValue - i*GRID_DIST

                    dist = math.sqrt((xMatrix-xCoord)**2 + (yMatrix-yCoord)**2)
                    if (dist < 2.5):  # 5% of GRID_DISTANCE
                        return gradx[i, j], grady[i, j]

                    if(distmin[0] == 0):  # first entry
                        distmin[0] = dist
                        p[0] = (xMatrix, yMatrix)
                        gradValues[0] = (gradx[i, j], grady[i, j])
                    elif(dist < distmin[0]):  # getting minimum distance
                        distmin[2] = distmin[1]
                        p[2] = p[1]
                        gradValues[2] = gradValues[1]
                        distmin[1] = distmin[0]
                        p[1] = p[0]
                        gradValues[1] = gradValues[0]
                        distmin[0] = dist
                        p[0] = (xMatrix, yMatrix)
                        gradValues[0] = gradx[i, j], grady[i, j]
                    elif(distmin[1] == 0):  # second entry
                        distmin[1] = dist
                        p[1] = (xMatrix, yMatrix)
                        gradValues[1] = gradx[i, j], grady[i, j]
                    elif(dist < distmin[1]):  # getting minimum distance
                        distmin[2] = distmin[1]
                        p[2] = p[1]
                        gradValues[2] = gradValues[1]
                        distmin[1] = dist
                        p[1] = (xMatrix, yMatrix)
                        gradValues[1] = gradx[i, j], grady[i, j]
                    elif(distmin[2] == 0 or dist < distmin[2]):  # third entry
                        distmin[2] = dist
                        p[2] = (xMatrix, yMatrix)
                        gradValues[2] = gradx[i, j], grady[i, j]

            # Loop that calculates sum of square distances
            for i, dist in enumerate(distmin):
                # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                if(dist < 0.5):
                    estimatedGradientX = gradValues[i][0]
                    estimatedGradientY = gradValues[i][1]
                    closePoint = True
                    continue
                else:
                    denom += 1/dist**2

            # if stick point is not so close to point on DEM:
            # there are no close points
            if(closePoint == False):
                for i, dist in enumerate(distmin):
                    # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                    estimatedGradientX += 1/dist**2/denom*gradValues[i][0]
                    estimatedGradientY += 1/dist**2/denom*gradValues[i][1]

            return estimatedGradientX, estimatedGradientY


def checkClosestPointSpeed(point):
    with open('vel_interpolada_3.dat', 'r') as velFile:
        p = None
        min = GRID_DIST
        v = None
        for line in velFile.readlines():
            line = line.split()

            xCoord, yCoord = float(line[0]), float(line[1])
            dist = math.sqrt((xCoord-point[0])**2 + (yCoord-point[1])**2)
            if(dist < min):
                min = dist
                p = (xCoord, yCoord)
                v = float(line[len(line)-1])

        return v


def updateStickPosition(gradX, gradY):
    stickPositions, grad = [], []
    sumaVel = []

    # Collects firsts positions
    with open('gradient_in_stick.dat', 'r') as stickPosAndGradient:
        for line in stickPosAndGradient.readlines():
            line = line.split()
            sumaVel.append(0.0)
            stickPositions.append((float(line[0]), float(line[1])))
            grad.append((float(line[len(line)-2]), float(line[len(line)-1])))

    with open('vel_darek.dat', 'r') as velFile:
        with open('vel_darek_mensuales.dat', 'w') as velMensuales:

            vel = np.loadtxt(velFile)
            initialDay = int(vel[0][len(vel[0])-4])

            for k, _ in enumerate(vel):
                day = int(vel[k][len(vel[k])-4])
                for i, _ in enumerate(stickPositions):

                    sumaVel[i] += vel[k][i]

                if (day - initialDay == 30):  # 30 days
                    velMensuales.write(f"{'  '.join(map(str,sumaVel))} \n")
                    initialDay = day
                    for i, _ in enumerate(stickPositions):
                        sumaVel[i] = 0.0

    with open('vel_darek_mensuales.dat', 'r') as velMensuales:
        with open('updated_stick_positions.dat', 'w') as updatedStickPositions:

            for line in velMensuales.readlines():
                line = line.split('  ')

                for i, pos in enumerate(stickPositions):
                    gradXvector, gradYvector = getGradientInSpecificStick(
                        gradX, gradY, (pos[0], pos[1]))
                    cosine = gradXvector / \
                        (math.sqrt(gradXvector**2 + gradYvector**2))
                    sine = gradYvector / \
                        (math.sqrt(gradXvector**2 + gradYvector**2))

                    if (line[i] == 'nan'):
                        # 21 = number of months between January 2013 to August 2014
                        xCoord = pos[0] + \
                            checkClosestPointSpeed((pos[0], pos[1]))/21*cosine
                        # 21 = number of months between January 2013 to August 2014
                        yCoord = pos[1] + \
                            checkClosestPointSpeed((pos[0], pos[1]))/21*sine
                    else:
                        xCoord = pos[0] + float(line[i])*cosine
                        yCoord = pos[1] + float(line[i])*sine
                    # speed * cosine(angle), where cosine(angle) = gradXvector/hypotenone
                    stickPositions[i] = (xCoord, yCoord)

                updatedStickPositions.write(
                    f" {'  '.join(map(str,stickPositions))} \n")


def checkClosestPoint(x, y):
    with open('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat', 'r') as demGrid:
        minDist = math.sqrt((x-xStartingValue)**2 + (y-yStartingValue)**2)
        for row in demGrid.readlines():
            row = row.split()
            dist = math.sqrt((x-float(row[0]))**2 + (y-float(row[1]))**2)
            if(dist < minDist):
                minDist = dist
                newX, newY = row[0], row[1]
        return newX, newY


def roundAndSelect():

    # selects contour points from DEM given contour points datafile (which are out of DEM)
    with open('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat', 'r') as demGrid:
        with open('zero_contour_points_Dem.dat', 'r') as file:
            with open('zero_contour_points_Dem_50_points.dat', 'w') as fileWrite:

                for line in file.readlines():
                    line = line.split()
                    x, y = round(float(line[0]), 0), round(float(line[1]), 0)
                    newX, newY = checkClosestPoint(x, y)
                    fileWrite.write(f"{newX} {newY} \n")


def printGradient(filenameX, filenameY):
    with open(f'{filenameX}.dat', 'r') as xgrad:
        with open(f'{filenameY}.dat', 'r') as ygrad:
            with open(f'{filenameX}_2plot.dat', 'w') as gradx:
                with open(f'{filenameY}_2plot.dat', 'w') as grady:

                    y = yFinalValue
                    for line in xgrad.readlines():
                        line = line.split()
                        x = xStartingValue
                        for i, _ in enumerate(line):
                            value = line[i]
                            gradx.write(f"{x} {y} {value} \n")
                            x += GRID_DIST
                        y -= GRID_DIST

                    y = yFinalValue
                    for line in ygrad.readlines():
                        line = line.split()
                        x = xStartingValue
                        for i, _ in enumerate(line):
                            value = line[i]
                            grady.write(f"{x} {y} {value} \n")
                            x += GRID_DIST
                        y -= GRID_DIST


# checks if given row,column are valid for given matrix
def checkArrayOutOfBoundOrZero(matrix, i, j):
    try:
        a = matrix[i][j]
        if (a == 0.0):
            return True
    except IndexError:
        return True
    return False


def customGradient(INTERVAL, MID_INTERVAL):

    distance = INTERVAL*2*GRID_DIST
    with open('cuadricula_superficie_Hansbreen.dat', 'r') as zGrid:
        with open(f'custom_gradient_x_{distance}m_interval.dat', 'w') as gradx:
            with open(f'custom_gradient_y_{distance}m_interval.dat', 'w') as grady:

                zMatrix = np.loadtxt(zGrid)

                gradXcustom, gradYcustom = np.zeros(
                    (329, 173), dtype=float), np.zeros((329, 173), dtype=float)

                # INTERVAL = 2 # This represents 2*GRID_DIST = 100m
                # MID_INTERVAL = int(INTERVAL/2)
                # For Y values
                for i, _ in enumerate(zMatrix):
                    for j, _ in enumerate(zMatrix[i]):

                        if(checkArrayOutOfBoundOrZero(zMatrix, i, j)):
                            gradYcustom[i][j] = 0.0
                            gradXcustom[i][j] = 0.0

                        # Y-gradient
                        # NOTE: - sign at the beginning because of grid coordinate system (from high y values to low y values)
                        if(checkArrayOutOfBoundOrZero(zMatrix, i-INTERVAL, j)):  # row above is empty or zero
                            # row immediately above is empty or zero
                            if(checkArrayOutOfBoundOrZero(zMatrix, i-MID_INTERVAL, j)):
                                if(checkArrayOutOfBoundOrZero(zMatrix, i+INTERVAL, j)):
                                    # next row empty or zero
                                    if(checkArrayOutOfBoundOrZero(zMatrix, i+MID_INTERVAL, j)):
                                        # no data because every relevant row is empty or zero
                                        gradYcustom[i][j] = 0.0
                                    else:
                                        gradYcustom[i][j] = -(zMatrix[i][j] - zMatrix[i+MID_INTERVAL][j])/(
                                            MID_INTERVAL*GRID_DIST)
                                else:
                                    gradYcustom[i][j] = -(zMatrix[i][j] -
                                                          zMatrix[i+INTERVAL][j])/(INTERVAL*GRID_DIST)

                            else:
                                if(checkArrayOutOfBoundOrZero(zMatrix, i+INTERVAL, j)):
                                    # next row empty or zero
                                    if(checkArrayOutOfBoundOrZero(zMatrix, i+MID_INTERVAL, j)):
                                        gradYcustom[i][j] = -(zMatrix[i-MID_INTERVAL][j] - zMatrix[i][j])/(
                                            MID_INTERVAL*GRID_DIST)
                                    else:
                                        gradYcustom[i][j] = -(zMatrix[i-MID_INTERVAL][j] - zMatrix[i+MID_INTERVAL][j])/(
                                            (MID_INTERVAL+MID_INTERVAL)*GRID_DIST)
                                else:
                                    gradYcustom[i][j] = -(zMatrix[i-MID_INTERVAL][j] - zMatrix[i+INTERVAL][j])/(
                                        (MID_INTERVAL+INTERVAL)*GRID_DIST)
                        else:
                            if(checkArrayOutOfBoundOrZero(zMatrix, i+INTERVAL, j)):
                                # next row empty or zero
                                if(checkArrayOutOfBoundOrZero(zMatrix, i+MID_INTERVAL, j)):
                                    gradYcustom[i][j] = -(zMatrix[i-INTERVAL]
                                                          [j] - zMatrix[i][j])/(INTERVAL*GRID_DIST)
                                else:
                                    gradYcustom[i][j] = -(zMatrix[i-INTERVAL][j] - zMatrix[i+MID_INTERVAL][j])/(
                                        (INTERVAL+MID_INTERVAL)*GRID_DIST)
                            else:
                                # general expression for gradient when all values are correct
                                gradYcustom[i][j] = -(zMatrix[i-INTERVAL][j] - zMatrix[i+INTERVAL][j])/(
                                    (INTERVAL+INTERVAL)*GRID_DIST)

                        #----------------------------------------------#
                        # X-gradient
                        # column left is empty or zero
                        if(checkArrayOutOfBoundOrZero(zMatrix, i, j-INTERVAL)):
                            # column immediately left is empty or zero
                            if(checkArrayOutOfBoundOrZero(zMatrix, i, j-MID_INTERVAL)):
                                if(checkArrayOutOfBoundOrZero(zMatrix, i, j+INTERVAL)):
                                    # next column empty or zero
                                    if(checkArrayOutOfBoundOrZero(zMatrix, i, j+MID_INTERVAL)):
                                        # no data because every relevant column is empty or zero
                                        gradXcustom[i][j] = 0.0
                                    else:
                                        gradXcustom[i][j] = (
                                            zMatrix[i][j] - zMatrix[i][j+MID_INTERVAL])/(MID_INTERVAL*GRID_DIST)
                                else:
                                    gradXcustom[i][j] = (
                                        zMatrix[i][j] - zMatrix[i][j+INTERVAL])/(INTERVAL*GRID_DIST)

                            else:
                                if(checkArrayOutOfBoundOrZero(zMatrix, i, j+INTERVAL)):
                                    # next column empty or zero
                                    if(checkArrayOutOfBoundOrZero(zMatrix, i, j+MID_INTERVAL)):
                                        gradXcustom[i][j] = (
                                            zMatrix[i][j-MID_INTERVAL] - zMatrix[i][j])/(MID_INTERVAL*GRID_DIST)
                                    else:
                                        gradXcustom[i][j] = (zMatrix[i][j-MID_INTERVAL] - zMatrix[i][j+MID_INTERVAL])/(
                                            (MID_INTERVAL+MID_INTERVAL)*GRID_DIST)
                                else:
                                    gradXcustom[i][j] = (
                                        zMatrix[i][j-MID_INTERVAL] - zMatrix[i][j+INTERVAL])/((MID_INTERVAL+INTERVAL)*GRID_DIST)
                        else:
                            if(checkArrayOutOfBoundOrZero(zMatrix, i, j+INTERVAL)):
                                # next column empty or zero
                                if(checkArrayOutOfBoundOrZero(zMatrix, i, j+MID_INTERVAL)):
                                    gradXcustom[i][j] = (
                                        zMatrix[i][j-INTERVAL] - zMatrix[i][j])/(INTERVAL*GRID_DIST)
                                else:
                                    gradXcustom[i][j] = (
                                        zMatrix[i][j-INTERVAL] - zMatrix[i][j+MID_INTERVAL])/((INTERVAL+MID_INTERVAL)*GRID_DIST)
                            else:
                                # general expression for gradient when all values are correct
                                gradXcustom[i][j] = (
                                    zMatrix[i][j-INTERVAL] - zMatrix[i][j+INTERVAL])/((INTERVAL+INTERVAL)*GRID_DIST)

                np.savetxt(gradx, gradXcustom, fmt='%.4f')
                np.savetxt(grady, gradYcustom, fmt='%.4f')


def gradientSofter(gradY, softMetrics=0):
    with open(f"{gradY}.dat", 'r') as grady2soft:
        with open(f"{gradY}_soft.dat", 'w') as gradySofted:
            grady = np.loadtxt(grady2soft)
            gradYcustom = np.zeros((329, 173), dtype=float)

            for i, _ in enumerate(grady):
                for j, value in enumerate(grady[i]):
                    sumGrad = value
                    # if (value > 0):
                    boolCheckHorizontal = {
                        j-2: checkArrayOutOfBoundOrZero(grady, i, j-2),
                        j-1: checkArrayOutOfBoundOrZero(grady, i, j-1),
                        j+1: checkArrayOutOfBoundOrZero(grady, i, j+1),
                        j+2: checkArrayOutOfBoundOrZero(grady, i, j+2)
                    }
                    boolCheck = {
                        i-2: checkArrayOutOfBoundOrZero(grady, i-2, j),
                        i-1: checkArrayOutOfBoundOrZero(grady, i-1, j),
                        i+1: checkArrayOutOfBoundOrZero(grady, i+1, j),
                        i+2: checkArrayOutOfBoundOrZero(grady, i+2, j)
                    }
                    # must be average between proper value and (4) more, so it starts from 1
                    count = 1
                    for key, values in boolCheck.items():
                        if(values == False):
                            sumGrad += grady[key, j]
                            count += 1
                    for key, values in boolCheckHorizontal.items():
                        if(values == False):
                            sumGrad += grady[i, key]
                            count += 1
                    gradYcustom[i][j] = sumGrad/count

            np.savetxt(gradySofted, gradYcustom, fmt='%.4f')


def stickPositionsPlotter():
    with open('updated_stick_positions.dat', 'r') as stickPositions:
        with open('updated_stick_position_2plot.dat', 'w') as stickPositions2plot:
            for line in stickPositions.readlines():
                for i, item in enumerate(line):
                    if (i == 0):
                        stickPositions2plot.write("")

                    elif (item == "(" or item == ")"):
                        stickPositions2plot.write(",")
                    else:
                        stickPositions2plot.write(item)


def yearlyStickSpeed():

    speeds = []
    for i in range(0, 16):
        speeds.append(0.0)

    monthCount = 0
    with open('vel_darek_mensuales.dat', 'r') as monthlySpeed:
        with open('vel_darek_anuales.dat', 'w') as yearlySpeed:
            for line in monthlySpeed.readlines():
                line = line.split('  ')
                monthCount += 1
                for i, _ in enumerate(speeds):
                    if (line[i] == 'nan'):
                        speeds[i] = -0.0
                    else:
                        speeds[i] = speeds[i] + float(line[i])
                if(monthCount == 12):
                    yearlySpeed.write(f" {'  '.join(map(str,speeds))} \n")
                    monthCount = 0


def monthlySticksResidue():
    yearAvgSpeeds, residues = np.zeros((6,16),dtype=float), np.zeros((73,16), dtype=float)
    
    with open('vel_darek_anuales.dat', 'r') as yearlySpeed:
        with open('vel_darek_mensuales.dat', 'r') as monthlySpeed:
            with open('residuos_mensuales.dat', 'w') as monthlyResidues:
                yearCount = 0
                for row in yearlySpeed.readlines():
                    row = row.split('  ')
                    for j, _ in enumerate(yearAvgSpeeds[0]):
                        if (row[j] == 'nan' or sign(float(row[j]))==-1.0):
                            yearAvgSpeeds[yearCount,j] = np.nan
                        else:
                            
                            yearAvgSpeeds[yearCount,j] = float(row[j])/12
                    yearCount += 1 
                
                month = 0
                for line in monthlySpeed.readlines():
                    line = line.split('  ')
                    year = 0
                    for j, _ in enumerate(residues[0]):
                        if((month+1) % 12 == 0):
                            year = (month+1)//12-1
                            if (line[j] == 'nan'):
                                residues[month,j] = np.nan
                            else:
                                residues[month,j] = float(line[j])-yearAvgSpeeds[year,j]
                        elif (line[j] == 'nan'):
                            residues[month,j] = np.nan
                        else:
                            residues[month,j] = float(line[j])-yearAvgSpeeds[year,j]
                    print(month)
                    month += 1
            
                np.savetxt(monthlyResidues, residues, fmt='%.4f')

def componentSplitter():
    year = 2005
    for month in range(5,78):
        with open(f'speedsAfterKrig/speedsAfterBayesianKriging_{month}_{year}.dat', 'r') as speeds:
            with open('custom_gradient_x_600m_interval_2plot.dat', 'r') as xGrad:
                with open('custom_gradient_y_600m_interval_2plot.dat', 'r') as yGrad:
                    with open(f'speedsAfterKrig/components/speedXComponentAfterBayesianKriging_{month}_{year}_.dat', 'w') as xSpeedFile:
                        with open(f'speedsAfterKrig/components/speedYComponentAfterBayesianKriging_{month}_{year}.dat', 'w') as ySpeedFile: 
                    
                            speeds = np.loadtxt(speeds)
                            xGrad = np.loadtxt(xGrad)
                            yGrad = np.loadtxt(yGrad)
                            for i, _ in enumerate(yGrad):
                                hyp = math.sqrt(xGrad[i][2]**2+yGrad[i][2]**2)
                                if(hyp == 0.0):
                                    xSpeed = 0.0
                                    ySpeed = 0.0
                                else:
                                    cosine = xGrad[i][2]/hyp
                                    sine = yGrad[i][2]/hyp
                                    xSpeed = speeds[i][2]*cosine
                                    ySpeed = speeds[i][2]*sine
                                xSpeedFile.write(f"{speeds[i][0]} {speeds[i][1]} {xSpeed} \n")
                                ySpeedFile.write(f"{speeds[i][0]} {speeds[i][1]} {ySpeed} \n")
        print(f'MONTH={month} YEAR={year}')
        if(month % 12 == 0):
            year +=1

def main():

    metros_gradiente = 600
    interval = int(metros_gradiente/(2*GRID_DIST))
    # customGradient(interval, interval-1)
    gradXfile = f"custom_gradient_x_{metros_gradiente}m_interval"
    gradYfile = f"custom_gradient_y_{metros_gradiente}m_interval"

    gradXfileDat = f"custom_gradient_x_{metros_gradiente}m_interval.dat"
    gradYfileDat = f"custom_gradient_y_{metros_gradiente}m_interval_soft.dat"

    # gradient(): this is conventional gradient
    # speedComponentsDem(gradXfileDat, gradYfileDat)
    # gradientSofter(gradYfile)
    # printGradient(gradXfile, f"{gradYfile}_soft")
    # getGradientInStickCoordinates(gradXfileDat, gradYfileDat)
    # updateStickPosition(gradXfileDat, gradYfileDat)
    # stickPositionsPlotter()

    softMetrics = metros_gradiente/(2*GRID_DIST)
    # monthlySticksResidue()
    componentSplitter()
    # yearlyStickSpeed()
    # roundAndSelect()
    # interpolate()
    

if __name__ == '__main__':
    main()
