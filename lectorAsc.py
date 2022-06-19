import numpy as np
xStartingValue, xFinalValue = 510425.000, 519025.000
yStartingValue, yFinalValue = 8547775.000, 8564175.000
GRID_DIST = 50


with open('output_residuos_2008_mes1.asc', 'r') as kriggedResInDem:
    with open('output_residuos_2008_mes1_2plot.dat', 'w') as plottedRes:
        residuosDem=np.loadtxt(kriggedResInDem)
        currentY = yFinalValue
        for i, _ in enumerate(residuosDem):
            currentX = xStartingValue
            for j, res in enumerate(residuosDem[i]):
                plottedRes.write(f"{currentX} {currentY} {res} \n")
                currentX += GRID_DIST
            currentY -= GRID_DIST


