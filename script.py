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

def getSpeedMap():
    with open ('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat','r') as demGrid:
        with open ('vel_interpolada.dat','r') as velGrid:
            with open ('velocidades_superficie_Dem.dat', 'w') as sortedGrid:

                    for line in demGrid.readlines():
                        value= line.split()
                        yIndex = yStartingValue
                        y = int((float(value[0]) - xStartingValue)/GRID_DIST)
                        x = 328 - int((float(value[1]) - yStartingValue)/GRID_DIST)
                        zDemGrid[x, y] = value[2]
                    for row, value in enumerate(zDemGrid):
                        sortedGrid.write(f"{' '.join(map(str,value))}\n")
            

posEstaca = []
with open ('posiciones_tyczki.dat','r') as posicionesEstacas:
    for line in posicionesEstacas.readlines():
        posEstaca.append(line.split()[0])
        posEstaca.append(line.split()[1])
        posEstaca.append(line.split()[2])
        posEstaca.append("XXXXXXXXXX")


def gradient():
    #Check gradient output 
    gradx, grady = np.gradient(zDemGrid)
    with open("gradient_x_dem.dat",'w') as gradXFile:
        with open("gradient_y_dem.dat",'w') as gradYFile:
            np.savetxt(gradXFile,gradx, fmt='%.4f')
            np.savetxt(gradYFile,grady, fmt='%.4f')
        
    # print(f"{gradx} \n {grady}")

def angleCalc():
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

                           

                
    




distMatrix = {
    "p1": [0,(0,0)],
    "p2": [0,(0,0)],
    "p3": [0,(0,0)]
}


        
    
def main():
    
    # interpolate()
    # gradient()
    angleCalc()
    # getSpeedMap()

    # fig = plt.figure(figsize=(6,6))
    # ax = fig.add_subplot(111, projection='3d')
    
    # # Plot a 3D surface
    # ax.plot_surface(np.array(xDemGrid), np.array(yDemGrid), Z )

    # plt.show()

    # interpolate()


if __name__=='__main__':
    main()
