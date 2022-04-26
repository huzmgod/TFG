import numpy as np
import math

##########################################################################################
###############################       CONSTANTS       ####################################
########################################################################################## 

###### STARTING GRID VALUES (HANSBREEN) WGS84 ######
xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue =  8547775.000, 8564175.000

###### GRID SQUARE LENGTH ######
GRID_DIST = 50

###### SIGN FUNCTION ######
sign = lambda x: math.copysign(1, x)

##########################################################################################
###############################         GRIDS         ####################################
########################################################################################## 

# Initial x-axis value: 510425
# Final x-axis value: 519025
# Initial y-axis value: 8547775
# Final y-axis value: 8564175

xDemGrid= np.linspace(510425.000,519025.000,num= 173)
yDemGrid= np.linspace(8547775.000, 8564175.000, num= 329)
zDemGrid =np.zeros((329,173), dtype=float)
xSpeedGrid, ySpeedGrid = [],[]
absSpeed = [] #m/year

#Data Storage: save speed values and dem grid from files 
with open ('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat','r') as demGrid:
    with open ('vel_recortada.dat','r') as velGrid:
            with open('cuadricula_superficie_Hansbreen.dat', 'w') as sortedGrid:
                
                for line in demGrid.readlines():
                    value= line.split()
                    yIndex = yStartingValue
                    y = int((float(value[0]) - xStartingValue)/GRID_DIST)
                    x = 328 - int((float(value[1]) - yStartingValue)/GRID_DIST)  # descendent order from greater y to lower y
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
    #Weighted average distance interpolation
    with open('vel_interpolada_3.dat','w') as outputSpeed:
        for i, x in enumerate(xDemGrid):
            print("change in X")
            for _, y in enumerate(yDemGrid):
                coordXdem=float(x)
                coordYdem=float(y)
                distmin=[0.0,0.0,0.0]
                speed=[0.0,0.0,0.0]
                p=[(0.0,0.0),(0.0,0.0),(0.0,0.0)] # vector of closest points
                denom = 0.0
                estimatedSpeed = 0.0
                closePoint = False
                print("change in Y")
                
                for k, _ in enumerate(xSpeedGrid):
                    coordXvel=float(xSpeedGrid[k])
                    coordYvel=float(ySpeedGrid[k])
                    modSpeed=float(absSpeed[k])
                    dist = math.sqrt((coordXdem-coordXvel)**2+(coordYdem-coordYvel)**2)
                    if(distmin[0]==0): #first entry
                        distmin[0]=dist
                        p[0]=(coordXvel,coordYvel)
                        speed[0]=modSpeed
                    elif(dist<distmin[0]): #getting minimum distance
                        distmin[2]=distmin[1]
                        p[2]=p[1]
                        speed[2]=speed[1]
                        distmin[1]=distmin[0]
                        p[1]=p[0]
                        speed[1]=speed[0]
                        distmin[0]=dist
                        p[0]=(coordXvel,coordYvel)
                        speed[0]=modSpeed
                    elif(distmin[1]==0): #second entry
                        distmin[1]=dist
                        p[1]=(coordXvel,coordYvel)
                        speed[1]=modSpeed
                    elif(dist<distmin[1]): #getting minimum distance
                        distmin[2]=distmin[1]
                        p[2]=p[1]
                        speed[2]=speed[1]
                        distmin[1]=dist
                        p[1]=(coordXvel,coordYvel)
                        speed[1]=modSpeed
                    elif(distmin[2]==0 or dist<distmin[2]): #third entry
                        distmin[2]=dist
                        p[2]=(coordXvel,coordYvel)
                        speed[2]=modSpeed

                              
                # Loop that calculates sum of square distances 
                for i, dist in enumerate(distmin):
                    # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                    if(dist <0.5):
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
                        # print(f"{x} {y} {distmin} {p}")
                
                outputSpeed.write(f"{coordXdem}  {coordYdem}  {p}  {distmin}  {speed}  {estimatedSpeed} \n")

def gradient():
    #Check gradient output 
    gradx, grady = np.gradient(zDemGrid)
    with open("gradient_x_dem.dat",'w') as gradXFile:
        with open("gradient_y_dem.dat",'w') as gradYFile:
            np.savetxt(gradXFile,gradx, fmt='%.4f')
            np.savetxt(gradYFile,grady, fmt='%.4f')
        
    # print(f"{gradx} \n {grady}")

def speedComponentsDem():
    #Calculate angle between gradient components
    with open ("vel_interpolada_3.dat",'r') as speedValues:
        with open("gradient_x_dem.dat",'r') as gradXFile:
            with open("gradient_y_dem.dat",'r') as gradYFile:
                    with open("speed_components.dat",'w') as speedFile:
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
                                        xValue = gradx[i,j]
                                        yValue = grady[i,j]
                                        xSpeed = 0.0
                                        ySpeed = 0.0
                                        speed = float(line.split('  ')[5]) 
                                        
                                        hypotenuse = math.sqrt(xValue**2 + yValue**2)
                                        
                                                                   
                                        if(xValue == 0 and yValue != 0):
                                            ySpeed = speed*sign(yValue)

                                        elif(yValue == 0 and xValue != 0):
                                            xSpeed = speed*sign(xValue)

                                        elif(xValue != 0 and yValue != 0):
                                            xSpeed = speed*xValue/hypotenuse
                                            ySpeed = speed*yValue/hypotenuse

                                        speedFile.write(f"{xCoord}  {yCoord}  {xMatrix}  {yMatrix}  {xSpeed}  {ySpeed}\n")
                                        break
                                if(equalCoords == True):
                                    break



def getGradientInStickCoordinates():

    ###### STICKS COORDINATES ######
    
    with open ('posiciones_tyczki.dat','r') as posicionesEstacas:
        with open ('gradient_x_dem.dat','r') as gradXFile:
            with open ('gradient_y_dem.dat', 'r') as gradYFile:
                with open ('gradient_in_stick.dat', 'w') as gradStickFile:
                    gradx = np.loadtxt(gradXFile)
                    grady = np.loadtxt(gradYFile)
                    
                    for line in posicionesEstacas.readlines():
                        distmin = [0.0,0.0,0.0]
                        xCoord = float(line.split()[0])
                        yCoord = float(line.split()[1])
                        p=[(0.0,0.0),(0.0,0.0),(0.0,0.0)] # vector of closest points
                        denom = 0.0
                        gradValues = [(0.0,0.0),(0.0,0.0),(0.0,0.0)]
                        closePoint = False
                        estimatedGradientX, estimatedGradientY = 0.0, 0.0
                        for i, _ in enumerate(gradx):
                            for j, _ in enumerate(gradx[i]):
                                xMatrix = xStartingValue + j*GRID_DIST
                                yMatrix = yFinalValue - i*GRID_DIST
                                dist = math.sqrt((xMatrix-xCoord)**2 + (yMatrix-yCoord)**2)     

                                if(distmin[0]==0): #first entry
                                    distmin[0]=dist
                                    p[0]=(xMatrix, yMatrix)
                                    gradValues[0]=(gradx[i,j],grady[i,j])
                                elif(dist<distmin[0]): #getting minimum distance
                                    distmin[2]=distmin[1]
                                    p[2]=p[1]
                                    gradValues[2]=gradValues[1]
                                    distmin[1]=distmin[0]
                                    p[1]=p[0]
                                    gradValues[1]=gradValues[0]
                                    distmin[0]=dist
                                    p[0]=(xMatrix, yMatrix)
                                    gradValues[0]=gradx[i,j],grady[i,j]
                                elif(distmin[1]==0): #second entry
                                    distmin[1]=dist
                                    p[1]=(xMatrix, yMatrix)
                                    gradValues[1]=gradx[i,j],grady[i,j]
                                elif(dist<distmin[1]): #getting minimum distance
                                    distmin[2]=distmin[1]
                                    p[2]=p[1]
                                    gradValues[2]=gradValues[1]
                                    distmin[1]=dist
                                    p[1]=(xMatrix, yMatrix)
                                    gradValues[1]=gradx[i,j],grady[i,j]
                                elif(distmin[2]==0 or dist<distmin[2]): #third entry
                                    distmin[2]=dist
                                    p[2]=(xMatrix, yMatrix)
                                    gradValues[2]=gradx[i,j],grady[i,j]

                                        
                        # Loop that calculates sum of square distances 
                        for i, dist in enumerate(distmin):
                            # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                            if(dist <0.5):
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
                                
                        
                        gradStickFile.write(f"{xCoord}  {yCoord}  {p}  {distmin}  {estimatedGradientX}  {estimatedGradientY} \n")


def getGradientInStickCoordinates(stickPosition):

    ###### This method calculates weighted average gradient around stick position ######
    
    with open ('gradient_x_dem.dat','r') as gradXFile:
        with open ('gradient_y_dem.dat', 'r') as gradYFile:
           
            gradx = np.loadtxt(gradXFile)
            grady = np.loadtxt(gradYFile)
        
            distmin = [0.0,0.0,0.0]

            xCoord = stickPosition[0]
            yCoord = stickPosition[1]

            p=[(0.0,0.0),(0.0,0.0),(0.0,0.0)] # vector of closest points
            denom = 0.0
            gradValues = [(0.0,0.0),(0.0,0.0),(0.0,0.0)]
            closePoint = False
            estimatedGradientX, estimatedGradientY = 0.0, 0.0

            for i, _ in enumerate(gradx):
                for j, _ in enumerate(gradx[i]):
                    xMatrix = xStartingValue + j*GRID_DIST
                    yMatrix = yFinalValue - i*GRID_DIST
                    dist = math.sqrt((xMatrix-xCoord)**2 + (yMatrix-yCoord)**2)     

                    if(distmin[0]==0): #first entry
                        distmin[0]=dist
                        p[0]=(xMatrix, yMatrix)
                        gradValues[0]=(gradx[i,j],grady[i,j])
                    elif(dist<distmin[0]): #getting minimum distance
                        distmin[2]=distmin[1]
                        p[2]=p[1]
                        gradValues[2]=gradValues[1]
                        distmin[1]=distmin[0]
                        p[1]=p[0]
                        gradValues[1]=gradValues[0]
                        distmin[0]=dist
                        p[0]=(xMatrix, yMatrix)
                        gradValues[0]=gradx[i,j],grady[i,j]
                    elif(distmin[1]==0): #second entry
                        distmin[1]=dist
                        p[1]=(xMatrix, yMatrix)
                        gradValues[1]=gradx[i,j],grady[i,j]
                    elif(dist<distmin[1]): #getting minimum distance
                        distmin[2]=distmin[1]
                        p[2]=p[1]
                        gradValues[2]=gradValues[1]
                        distmin[1]=dist
                        p[1]=(xMatrix, yMatrix)
                        gradValues[1]=gradx[i,j],grady[i,j]
                    elif(distmin[2]==0 or dist<distmin[2]): #third entry
                        distmin[2]=dist
                        p[2]=(xMatrix, yMatrix)
                        gradValues[2]=gradx[i,j],grady[i,j]

                            
            # Loop that calculates sum of square distances 
            for i, dist in enumerate(distmin):
                # if distmin[0] < 0.5 => dist ~= 0, then no need to interpolate the grid
                if(dist <0.5):
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
    with open ('vel_interpolada_3.dat', 'r') as velFile:
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




def updateStickPosition():
    stickPositions, grad = [], [], []
    sumaVel = []
    
    with open ('gradient_in_stick.dat', 'r') as stickPosAndGradient:
        for line in stickPosAndGradient.readlines():
            line = line.split()
            sumaVel.append(0.0)
            stickPositions.append((float(line[0]), float(line[1])))        
            grad.append((float(line[len(line)-2]),float(line[len(line)-1])))

    with open ('vel_darek.dat', 'r') as velFile:
        with open ('vel_darek_mensuales.dat', 'w') as velMensuales:
            
            vel = np.loadtxt(velFile)
            initialDay = int(vel[0][len(vel[0])-4])
            
            for k, _ in enumerate(vel):
                day = int(vel[k][len(vel[k])-4])
                for i, _ in enumerate(stickPositions):
                    

                    sumaVel[i] += vel[k][i]
                    # print(sumaVel[i])
                if (day - initialDay ==30): #30 days
                    velMensuales.write(f"{'  '.join(map(str,sumaVel))} \n")
                    initialDay = day
                    for i, _ in enumerate(stickPositions):
                        sumaVel[i] = 0.0

    with open ('vel_darek_mensuales.dat', 'r') as velMensuales:
        with open ('updated_stick_positions.dat', 'w') as updatedStickPositions:
            
            
            for line in velMensuales.readlines():
                line = line.split('  ')
                
                for i, pos in enumerate(stickPositions):
                    gradXvector, gradYvector = getGradientInStickCoordinates((pos[0],pos[1]))
                    
                    cosine = gradXvector/(math.sqrt(gradXvector**2 + gradYvector**2))
                    sine = gradYvector/(math.sqrt(gradXvector**2 + gradYvector**2))

                    if (line[i] == 'nan'):
                        print("nan")
                        xCoord = pos[0] + checkClosestPointSpeed((pos[0],pos[1]))/12*cosine
                        yCoord = pos[1] + checkClosestPointSpeed((pos[0],pos[1]))/12*sine
                    else:
                        xCoord = pos[0] + float(line[i])*cosine
                        yCoord = pos[1] +float(line[i])*sine
                    stickPositions[i] = (xCoord,yCoord)  # speed * cosine(angle), where cosine(angle) = gradXvector/hypotenone 
                    
                
                updatedStickPositions.write(f" {'  '.join(map(str,stickPositions))} \n")





def main():
    
    # interpolate()
    # gradient()
    # speedComponentsDem()
    # getSpeedMap()
    # getGradientInStickCoordinates()
    updateStickPosition()
    # fig = plt.figure(figsize=(6,6))
    # ax = fig.add_subplot(111, projection='3d')
    
    # # Plot a 3D surface
    # ax.plot_surface(np.array(xDemGrid), np.array(yDemGrid), Z )

    # plt.show()

    # interpolate()


if __name__=='__main__':
    main()
