#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:54:53 2023

@author: jaime
"""

import numpy as np
import math
import datetime
import ClsKriging_BK as kg

##########################################################################################
###############################       CONSTANTS       ####################################
##########################################################################################


###### TEMPORAL STEP #######
TIME_STEP = 30.4375

###### GRID SQUARE LENGTH ######
GRID_DIST = 50

###### STARTING GRID VALUES (HANSBREEN) WGS84 ######
xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
rows = (yFinalValue-yStartingValue)/GRID_DIST
cols = (xFinalValue-xStartingValue)/GRID_DIST

##### CALENDAR DICTIONARY #####
'''
Code below generates a calendar dictionary:
'''


def generateCalendarDict():
    start_month = 5
    start_year = 2005
    num_months = 72

    # Function to generate month name from a given number
    def month_name(month_number):
        return datetime.date(1900, month_number, 1).strftime('%B').lower()

    # Function to generate dictionary keys and values
    def generate_month_year_dict(start_month, start_year, num_months):
        month_year_dict = {}

        current_month = start_month
        current_year = start_year

        for i in range(1, num_months + 1):
            key = f"{current_year}_{current_month}"
            month_year_dict[str(i)] = key

            current_month += 1
            if current_month > 12:
                current_month = 1
                current_year += 1

        return month_year_dict

    # Generate the dictionary
    month_year_dict = generate_month_year_dict(
        start_month, start_year, num_months)
    return month_year_dict


calendar = {'1': '2005_5',
            '2': '2005_6',
            '3': '2005_7',
            '4': '2005_8',
            '5': '2005_9',
            '6': '2005_10',
            '7': '2005_11',
            '8': '2005_12',
            '9': '2006_1',
            '10': '2006_2',
            '11': '2006_3',
            '12': '2006_4',
            '13': '2006_5',
            '14': '2006_6',
            '15': '2006_7',
            '16': '2006_8',
            '17': '2006_9',
            '18': '2006_10',
            '19': '2006_11',
            '20': '2006_12',
            '21': '2007_1',
            '22': '2007_2',
            '23': '2007_3',
            '24': '2007_4',
            '25': '2007_5',
            '26': '2007_6',
            '27': '2007_7',
            '28': '2007_8',
            '29': '2007_9',
            '30': '2007_10',
            '31': '2007_11',
            '32': '2007_12',
            '33': '2008_1',
            '34': '2008_2',
            '35': '2008_3',
            '36': '2008_4',
            '37': '2008_5',
            '38': '2008_6',
            '39': '2008_7',
            '40': '2008_8',
            '41': '2008_9',
            '42': '2008_10',
            '43': '2008_11',
            '44': '2008_12',
            '45': '2009_1',
            '46': '2009_2',
            '47': '2009_3',
            '48': '2009_4',
            '49': '2009_5',
            '50': '2009_6',
            '51': '2009_7',
            '52': '2009_8',
            '53': '2009_9',
            '54': '2009_10',
            '55': '2009_11',
            '56': '2009_12',
            '57': '2010_1',
            '58': '2010_2',
            '59': '2010_3',
            '60': '2010_4',
            '61': '2010_5',
            '62': '2010_6',
            '63': '2010_7',
            '64': '2010_8',
            '65': '2010_9',
            '66': '2010_10',
            '67': '2010_11',
            '68': '2010_12',
            '69': '2011_1',
            '70': '2011_2',
            '71': '2011_3',
            '72': '2011_4'}

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

def interpolateDEMSpeedsIntoPoint(x, y, dataPath="datafiles/vel_y_cambiada_blanked_no0.dat"):
    '''
    This function interpolates DEM speed values into a given point using IDW.
    x: float, X coordinate of the given point
    y: float, Y coordinate of the given point
    dataPath: string, path to the file containing the DEM speed values. By default, satellite interpolated in DEM speed values.
    Returns monthly speed value (estimatedSpeed / 12)
    '''
    speedDataPath = dataPath
    dem_data = np.loadtxt(speedDataPath)
    xDemGrid, yDemGrid, demSpeeds = dem_data[:,
                                             0], dem_data[:, 1], dem_data[:, 2]

    distmin = [0.0, 0.0, 0.0]
    speed = [0.0, 0.0, 0.0]
    p = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]
    denom = 0.0
    estimatedSpeed = 0.0
    closePoint = False

    for i, (coordXdem, coordYdem, modSpeed) in enumerate(zip(xDemGrid, yDemGrid, demSpeeds)):
        dist = math.sqrt((x - coordXdem) ** 2 + (y - coordYdem) ** 2)

        if (distmin[0] == 0):  # first entry
            distmin[0] = dist
            p[0] = (coordXdem, coordYdem)
            speed[0] = modSpeed
        elif (dist < distmin[0]):  # getting minimum distance
            distmin[2] = distmin[1]
            p[2] = p[1]
            speed[2] = speed[1]
            distmin[1] = distmin[0]
            p[1] = p[0]
            speed[1] = speed[0]
            distmin[0] = dist
            p[0] = (coordXdem, coordYdem)
            speed[0] = modSpeed
        elif (distmin[1] == 0):  # second entry
            distmin[1] = dist
            p[1] = (coordXdem, coordYdem)
            speed[1] = modSpeed
        elif (dist < distmin[1]):  # getting minimum distance
            distmin[2] = distmin[1]
            p[2] = p[1]
            speed[2] = speed[1]
            distmin[1] = dist
            p[1] = (coordXdem, coordYdem)
            speed[1] = modSpeed
        elif (distmin[2] == 0 or dist < distmin[2]):  # third entry
            distmin[2] = dist
            p[2] = (coordXdem, coordYdem)
            speed[2] = modSpeed

    for i, dist in enumerate(distmin):
        if (dist < 0.5):
            estimatedSpeed = speed[i]
            closePoint = True
            break
        else:
            denom += 1 / dist ** 2

    if (closePoint == False):
        for i, dist in enumerate(distmin):
            estimatedSpeed += 1 / dist ** 2 / denom * speed[i]

    return estimatedSpeed

def getGradientInSpecificStick(gradX, gradY, stickPosition):
    '''
    This method calculates weighted average gradient of stick position interpolating gradients of closest points
    '''

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

                    if (distmin[0] == 0):  # first entry
                        distmin[0] = dist
                        p[0] = (xMatrix, yMatrix)
                        gradValues[0] = (gradx[i, j], grady[i, j])
                    elif (dist < distmin[0]):  # getting minimum distance
                        distmin[2] = distmin[1]
                        p[2] = p[1]
                        gradValues[2] = gradValues[1]
                        distmin[1] = distmin[0]
                        p[1] = p[0]
                        gradValues[1] = gradValues[0]
                        distmin[0] = dist
                        p[0] = (xMatrix, yMatrix)
                        gradValues[0] = gradx[i, j], grady[i, j]
                    elif (distmin[1] == 0):  # second entry
                        distmin[1] = dist
                        p[1] = (xMatrix, yMatrix)
                        gradValues[1] = gradx[i, j], grady[i, j]
                    elif (dist < distmin[1]):  # getting minimum distance
                        distmin[2] = distmin[1]
                        p[2] = p[1]
                        gradValues[2] = gradValues[1]
                        distmin[1] = dist
                        p[1] = (xMatrix, yMatrix)
                        gradValues[1] = gradx[i, j], grady[i, j]
                    elif (distmin[2] == 0 or dist < distmin[2]):  # third entry
                        distmin[2] = dist
                        p[2] = (xMatrix, yMatrix)
                        gradValues[2] = gradx[i, j], grady[i, j]

            # Loop that calculates sum of square distances
            for i, dist in enumerate(distmin):
                # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                if (dist < 0.5):
                    estimatedGradientX = gradValues[i][0]
                    estimatedGradientY = gradValues[i][1]
                    closePoint = True
                    continue
                else:
                    denom += 1/dist**2

            # if stick point is not so close to point on DEM:
            # there are no close points
            if (closePoint == False):
                for i, dist in enumerate(distmin):
                    # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                    estimatedGradientX += 1/dist**2/denom*gradValues[i][0]
                    estimatedGradientY += 1/dist**2/denom*gradValues[i][1]

            return estimatedGradientX, estimatedGradientY

def fixedResidue2Kriging():
    '''
    For each month, calculates the residues between the monthly speed and the prior monthly speed (only sticks).
    Then, it writes the residue in a file to be used as input for ordinaryKrigingResidues.py, which interpolates the residues
    into the whole DEM.
    Right after, we add the interpolated residues to the prior speed values to get the final speed values.
    '''
    monthlySpeedValues = np.loadtxt('datafiles/vel_darek_mensuales_30.4375_normalized.dat')
    priorSpeedValues = np.loadtxt('datafiles/vel_interpolada_DEM_ordenada.dat')
    speedValues = np.zeros(priorSpeedValues.shape)

    residues = np.zeros(monthlySpeedValues.shape)
    # First positions (from data)
    positions = np.loadtxt('datafiles/posiciones_tyczki.dat')

    monthKey = 1

    def calculateResidues(monthKey, positions):
        # Case month 1
        if (monthKey == 1):
            for j, value in enumerate(monthlySpeedValues[monthKey-1]):
                if (np.isnan(value)):
                    residues[monthKey-1][j] = "NaN"
                else:
                    residues[monthKey-1][j] = value - \
                        interpolateDEMSpeedsIntoPoint(
                            positions[j][0], positions[j][1])/12

        # Case other months
        else:
            for j, value in enumerate(monthlySpeedValues[monthKey-1]):
                if (np.isnan(value)):
                    residues[monthKey-1][j] = "NaN"
                else:
                    residues[monthKey-1][j] = value - interpolateDEMSpeedsIntoPoint(positions[j][0], positions[j][1], \
                     f"speedsAfterKrig/fixedSpeedModules2/output_velocidades_{calendar[str(monthKey-1)]}_2plot.dat")

    def ordinaryKriging(positions, residues, monthKey):

        ##################   ORDINARY KRIGING (RESIDUES)   ###########################

        with open('datafiles/zero_contour_points_Dem_40_points.dat', 'r') as contourZeros:
            zeros = np.loadtxt(contourZeros)
            resMatrix = residues
            plot = False
            nlags = 56
            krigData, aux = [], []
            for i, _ in enumerate(zeros):
                inputData = [zeros[i][0], zeros[i][1], 0.0]
                krigData.append(inputData)
            for i in range(0, 16):
                coords = positions[i][0], positions[i][1]
                xCoord = float(coords[0])
                yCoord = float(coords[1])
                valueRes = resMatrix[monthKey-1][i]
                if (np.isnan(valueRes)
                        or np.isnan(xCoord)):  # associated to undefined stick speed
                    continue
                data = [xCoord, yCoord, valueRes]
                # print(f'DATA = {data}')
                krigData.append(data)
            krigData = np.array(krigData)

            dataPath="datafiles/vel_y_cambiada_blanked_no0.dat"

            dem_data = np.loadtxt(dataPath)
            xDemGrid, yDemGrid, demSpeeds = dem_data[:,0], dem_data[:, 1], dem_data[:, 2]
            DemGrid = np.zeros((len(xDemGrid),2))
            for i in range(len(xDemGrid)):
                DemGrid[i,0] = dem_data[i,0]
                DemGrid[i,1] = dem_data[i,1]
            #DemGrid = np.array(DemGrid)

            interpolate = kg.Kriging().ordinary(krigData, DemGrid)

            np.savetxt(f'speedsAfterKrig/fixedResidues2/output_residuos_{calendar[str(monthKey)]}_2plot.dat', interpolate) #, fmt=['%.1f', '%.1f', '%.15f'])

    def crossValidation(positions, residues, monthKey):

        ##################   ORDINARY KRIGING (RESIDUES)   ###########################

        with open('datafiles/zero_contour_points_Dem_40_points.dat', 'r') as contourZeros:
            zeros = np.loadtxt(contourZeros)
            resMatrix = residues
            nlags = 56
            krigData, aux = [], []
            for i, _ in enumerate(zeros):
                inputData = [zeros[i][0], zeros[i][1], 0.0]
                krigData.append(inputData)

            residuals = []
            for i in range(0, 16):
                currentKrigData = krigData.copy()

                coords = positions[i][0], positions[i][1]
                xCoord = float(coords[0])
                yCoord = float(coords[1])
                valueRes = resMatrix[monthKey-1][i]
                if (np.isnan(valueRes) or np.isnan(xCoord)):
                    continue

                for j in range(0, 16):
                    if j != i:
                        coords = positions[j][0], positions[j][1]
                        xCoord = float(coords[0])
                        yCoord = float(coords[1])
                        valueRes = resMatrix[monthKey-1][j]
                        if (np.isnan(valueRes) or np.isnan(xCoord)):
                            continue
                        data = [xCoord, yCoord, valueRes]
                        currentKrigData.append(data)
                currentKrigData = np.array(currentKrigData)

                #TODO: check how [xCoord, yCoord] is passed to Kriging()
                interpolate = kg.Kriging().ordinary(currentKrigData, [xCoord, yCoord])

                residual = valueRes - interpolate[2]
                residuals.append(residual)
            with open (f'crossValidationResidues', 'w') as outputFile:
                for i, _ in enumerate(residuals):
                    outputFile.write(f'{residuals[i]}\n')
            #return np.array(residuals)

    def bayesianKriging(monthKey):
        ######## BAYESIAN KRIGING (RESIDUES) ########
        with open(f'speedsAfterKrig/fixedResidues2/output_residuos_{calendar[str(monthKey)]}_2plot.dat', 'r') as kriggedResidues:
            kriggedResidues = np.loadtxt(kriggedResidues)
            if (monthKey == 1):
                for i, row in enumerate(speedValues):
                    speedValues[i][0] = priorSpeedValues[i][0]
                    speedValues[i][1] = priorSpeedValues[i][1]
                    speedValues[i][2] = priorSpeedValues[i][2]/12 + kriggedResidues[i][2]
                np.savetxt(
                    f'speedsAfterKrig/fixedSpeedModules2/output_velocidades_{calendar[str(monthKey)]}_2plot.dat', speedValues, fmt=['%.1f', '%.1f', '%.15f'])
            else:
                priorValues = np.loadtxt(f"speedsAfterKrig/fixedSpeedModules2/output_velocidades_{calendar[str(monthKey-1)]}_2plot.dat")
                for i, row in enumerate(speedValues):
                    speedValues[i][0] = priorValues[i][0]
                    speedValues[i][1] = priorValues[i][1]
                    speedValues[i][2] = priorValues[i][2] + kriggedResidues[i][2]
                np.savetxt(
                    f'speedsAfterKrig/fixedSpeedModules2/output_velocidades_{calendar[str(monthKey)]}_2plot.dat', speedValues, fmt=['%.1f', '%.1f', '%.15f'])

    # What have we done so far?
    # We have calculated the residues of the first month (May 2005) and we have krigged them into the DEM.
    # With those krigged residues, we have calculated the speed modules for the first month using the prior speed modules (satellite speeds).
    # Now we have to calculate the new stick positions and repeat the process for the following 71 months.

    def calculateNewStickPositions(monthKey, positions):
        '''
        Calculates the stick positions for the month (monthkey+1) using the speed modules just calculated.
        To do so, we get gradients from the gradient files and calculate the new position as p1_x = p0_x + v * g_x/(|g|), p1_y = p0_y + g_y/(|g|) * v
        '''
        gradXFileName = f'gradients/custom_gradient_x_600m_interval.dat'
        gradYFileName = f'gradients/custom_gradient_y_600m_interval_soft.dat'
        speedValues = np.loadtxt(
            f'speedsAfterKrig/fixedSpeedModules2/output_velocidades_{calendar[str(monthKey)]}_2plot.dat')
        newPositions = np.zeros(positions.shape)
        for i, _ in enumerate(positions):
            coords = positions[i][0], positions[i][1]
            gradientValues = getGradientInSpecificStick(
                gradXFileName, gradYFileName, coords)
            newPositions[i][0] = positions[i][0] + speedValues[i][2] * \
                gradientValues[0]/(math.sqrt(gradientValues[0]
                                   ** 2 + gradientValues[1]**2))
            newPositions[i][1] = positions[i][1] + speedValues[i][2] * \
                gradientValues[1]/(math.sqrt(gradientValues[0]
                                   ** 2 + gradientValues[1]**2))
        return newPositions


    # Now we automate the process for the 72 months.
    storedPositions = [[(0, 0) for i in range(16)] for j in range(72)]
    for monthKey in range(1, 73):
        for i, row in enumerate(positions):
            x, y = row[0], row[1]
            storedPositions[monthKey-1][i] = x, y
        calculateResidues(monthKey, positions)
        ordinaryKriging(positions, residues, monthKey)
        bayesianKriging(monthKey)
        positions = calculateNewStickPositions(monthKey, positions)
        crossValidation(positions, residues, monthKey)

    """ with open('storedPositions.dat', 'w') as storedPositionsFile:
        for i, row in enumerate(storedPositions):
            for j, stick in enumerate(row):
                storedPositionsFile.write(
                    f"{storedPositions[i][j][0]} {storedPositions[i][j][1]}, ")
            storedPositionsFile.write("\n") """

def main():

    metros_gradiente = 600
    interval = int(metros_gradiente/(2*GRID_DIST))
    # customGradient(interval, interval-1)

    gradXfile = f"gradients/custom_gradient_x_{metros_gradiente}m_interval"
    gradYfile = f"gradients/custom_gradient_y_{metros_gradiente}m_interval"

    gradXfileDat = f"gradients/custom_gradient_x_{metros_gradiente}m_interval.dat"
    gradYfileDat = f"gradients/custom_gradient_y_{metros_gradiente}m_interval_soft.dat"

    softMetrics = metros_gradiente/(2*GRID_DIST)

    # speedComponentsDem(gradXfileDat, gradYfileDat)
    # gradientSofter(gradYfile)
    # printGradient(gradXfile, f"{gradYfile}_soft")
    # getGradientInStickCoordinates(gradXfileDat, gradYfileDat)
    # computeMonthlySpeeds()
    # updateStickPosition(gradXfileDat, gradYfileDat)
    # stickPositionsPlotter()
    # yearlyStickSpeed()

    # monthlySticksResidue()
    # componentSplitter()

    # roundAndSelect()
    # interpolate()

    fixedResidue2Kriging()
    #componentSplitter()


if __name__ == '__main__':
    main()
