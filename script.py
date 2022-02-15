import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import math

xDemGrid, yDemGrid = [],[]
xSpeedGrid, ySpeedGrid = [],[]
absSpeed = [] #m/year
GRID_DIM = 40

#Data Storage
with open ('dem_utm_wgs1_clipped_Hans_rgi60_50m_recortado_blanked.dat','r') as demGrid:
    with open ('vel_recortada.dat','r') as velGrid:
        for line in demGrid.readlines():
            xDemGrid.append(line.split()[0])
            yDemGrid.append(line.split()[1])
        for line in velGrid.readlines():
            xSpeedGrid.append(line.split()[0])
            ySpeedGrid.append(line.split()[1])
            absSpeed.append(line.split()[3])

# grid = griddata((xSpeedGrid,ySpeedGrid), absSpeed, (xDemGrid,yDemGrid), method= 'cubic')
# with open ('vel_interpolada.dat','w') as interpolGrid:
#     for elem in list(grid.T):
#         interpolGrid.write(f"{elem}\n")
distMatrix = {
    "p1": [0,(0,0)],
    "p2": [0,(0,0)],
    "p3": [0,(0,0)]
}
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
            if(dist>GRID_DIM):
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
            denom += 1/distmin[i]**2
        for i in range(0,len(distmin)):
            estimatedSpeed += (1/distmin[i]**2) * (1/denom) * speed[i]
        outputSpeed.write(f"{coordXdem}  {coordYdem}  {estimatedSpeed}\n")
        
    