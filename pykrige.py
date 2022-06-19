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


# data = np.array(
#     [
#         [510425.000, 8547775.000, 0.47],
#         [510475.000, 8547825.000, 0.56],
#         [510525.000, 8547875.000, 0.74],
#         [510575.000, 8547925.000, 1.47],
#         [510625.000, 8547975.000, 1.74],
#     ]
# )
krigData = []

with open('zero_contour_points_Dem_50_points_random.dat', 'r') as contourZeros:
    zeros = np.loadtxt(contourZeros)
    for i, _ in enumerate(zeros):
        data = [zeros[i][0],zeros[i][1],0.0]
        krigData.append(data)
 
with open('updated_stick_position_2plot.dat', 'r') as stickPositions:
    with open('residuos_mensuales.dat', 'r') as residues:
        res = np.loadtxt(residues)
        month = 0
        for position in stickPositions.readlines():
            aux = []
            positions = position.split(',  ,')
            for i in range(0, 16):
                stick1 = str(positions[i].split(', '))
                for char in ['[', ']', ' ', "'", "\n"]:
                    stick1 = stick1.replace(char, "")
                coords = stick1.split(',')
                xCoord = float(coords[0])
                yCoord = float(coords[1])
                valueRes = res[month, i]
                data = [xCoord, yCoord, valueRes]

                krigData.append(data)
            month += 1
            break
krigData = np.array(krigData)
print(krigData)
###############################################################################
# Create the ordinary kriging object. Required inputs are the X-coordinates of
# the data points, the Y-coordinates of the data points, and the Z-values of the
# data points. If no variogram model is specified, defaults to a linear variogram
# model. If no variogram model parameters are specified, then the code automatically
# calculates the parameters by fitting the variogram model to the binned
# experimental semivariogram. The verbose kwarg controls code talk-back, and
# the enable_plotting kwarg controls the display of the semivariogram.

OK = OrdinaryKriging(
    krigData[:, 0],
    krigData[:, 1],
    krigData[:, 2],
    variogram_model="gaussian",
    verbose=False,
    enable_plotting=False,
    enable_statistics=True,
    exact_values = True
)

###############################################################################
# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points.
# (See OrdinaryKriging.__doc__ for more information.)

z, ss = OK.execute("grid", xDemGrid, yDemGrid)

# ###############################################################################
# # Writes the kriged grid to an ASCII grid file and plot it.

kt.write_asc_grid(xDemGrid, yDemGrid, z, filename="output.asc")
plt.imshow(z)
plt.show()
