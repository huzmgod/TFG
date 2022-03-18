import numpy as np
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import math

xDemGrid, yDemGrid, zDemGrid = [],[],np.zeros((173,329))
xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
xSpeedGrid, ySpeedGrid = [],[]
absSpeed = [] #m/year
GRID_DIST = 50

##########################################################################################
###############################         GRIDS         ####################################
########################################################################################## 
# Initial x-axis value: 510425
# Final x-axis value: 519025
# Initial y-axis value: 8547775
# Final y-axis value: 8564175

xDem= np.linspace(510425.000,519025.000,num= 173)
yDem= np.linspace(8547775.000, 8564175.000, num= 329)

#Data Storage: save speed values and dem grid from files 
with open ('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat','r') as demGrid:
    with open ('vel_recortada.dat','r') as velGrid:
            with open('cuadricula_superficie_Hansbreen.dat', 'w') as sortedGrid:
                
                for line in demGrid.readlines():
                    value= line.split()
                    xIndex = xStartingValue
                    for row in range(len(zDemGrid)): #for every x value
                        yIndex = yStartingValue
                        if (float(value[0])==xIndex):
                            for column in range(len(zDemGrid[row])):    #for every y value
                                if (float(value[1])==yIndex):
                                    zDemGrid[row,column] = float(value[2])
                                yIndex += GRID_DIST
                        xIndex += GRID_DIST
                for row in zDemGrid:
                    np.savetxt(sortedGrid,row,fmt='%2f')
                    #Not needed anymore. These are the lists of possible values (grids x,y,z)
                            # if value[0] not in xDemGrid:
                            #     xDemGrid.append(value[0])
                            # if value[1] not in yDemGrid:
                            #     yDemGrid.append(value[1])
                            # zDemGrid= np.vstack([zDemGrid, line.split()[2]])
                            
                for line in velGrid.readlines():
                    xSpeedGrid.append(line.split()[0])
                    ySpeedGrid.append(line.split()[1])
                    absSpeed.append(line.split()[3])



posEstaca = []
with open ('posiciones_tyczki.dat','r') as posicionesEstacas:
    for line in posicionesEstacas.readlines():
        posEstaca.append((line.split()[0],line.split()[1],line.split()[2]))

def gradient():
    #Check gradient output 
    gradx, grady = np.gradient(zDemGrid, xDem, yDem)
    with open("gradient_dem.dat",'ab') as gradFile:
        np.savetxt(gradFile,gradx)
        gradFile.write(b"\n")
        np.savetxt(gradFile,grady)
        
    # print(f"{gradx} \n {grady}")


distMatrix = {
    "p1": [0,(0,0)],
    "p2": [0,(0,0)],
    "p3": [0,(0,0)]
}

def interpolate():
    #Weighted average distance interpolation
    with open('vel_interpolada.dat','w') as outputSpeed:
        for i, value in enumerate(xDemGrid):
            coordXdem=float(xDemGrid[i])
            coordYdem=float(yDemGrid[i])
            distmin=[0.0,0.0,0.0]
            speed=[0.0,0.0,0.0]
            p=[(0.0,0.0),(0.0,0.0),(0.0,0.0)]
            denom = 0.0
            estimatedSpeed = 0.0
            for j, _ in enumerate(xSpeedGrid):
                coordXvel=float(xSpeedGrid[j])
                coordYvel=float(ySpeedGrid[j])
                modSpeed=float(absSpeed[j])
                dist = math.sqrt((coordXdem-coordXvel)**2+(coordYdem-coordYvel)**2)
                if(dist>GRID_DIST):
                    continue
                elif(distmin[0]==0): #first entry
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
            for i in range(0,len(distmin)):
                # k is 2 in this algorithm
                if (distmin[i]**2-0<0.1):
                    denom += 1/distmin[i]
                else: denom += 1/distmin[i]**2
            for i in range(0,len(distmin)):
                if (distmin[i]==0):
                    estimatedSpeed = speed[i]
                    continue
                else:
                    estimatedSpeed += (1/distmin[i]**2) * (1/denom) * speed[i]
            outputSpeed.write(f"{coordXdem}  {coordYdem}  {estimatedSpeed}\n")
        
    
def main():
#    print(np.sort(xDemGrid).size)
#    print(np.sort(yDemGrid).size)
#    print(xDem)
#    print(yDem)
    gradient()
    # print(zDemGrid)
    print(posEstaca)

    # fig = plt.figure(figsize=(6,6))
    # ax = fig.add_subplot(111, projection='3d')
    
    # # Plot a 3D surface
    # ax.plot_surface(np.array(xDemGrid), np.array(yDemGrid), Z )

    # plt.show()

    # interpolate()


if __name__=='__main__':
    main()