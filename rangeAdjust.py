import random
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean
from scipy.optimize import minimize
import gstools as gs

# Load data
data = np.loadtxt("datafiles/vel_interpolada_DEM_ordenada.dat")

#Calculating the euclidean distance between the points
def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

#Extracting 1000 random points from the data
def extract_random_points(data, n):
    return random.sample(data, n)

#Calculating the experimental variogram
def experimental_variogram(points):
    variogram = []
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            x1, y1, v1 = points[i]
            x2, y2, v2 = points[j]
            distance = euclidean_distance(x1, y1, x2, y2)
            if distance >= 50 and distance <= 50000:
                variogram.append((distance, 0.5 * (v1 - v2)**2))
    variogram = sorted(variogram, key=lambda x: x[0])
    return variogram

#Fitting the experimental variogram to a spherical theoretical model
def spherical_model(distance, sill, range_, nugget=0):
    if distance <= range_:
        return nugget + (sill - nugget) * ((3 * distance)/(2 * range_) - (distance / range_)**3)
    else:
        return sill

def fit_spherical_model(variogram, sill, range_, nugget=0):
    fit = []
    for distance, gamma in variogram:
        fit.append((distance, spherical_model(distance, sill, range_, nugget)))
    return fit

#Repeating the process 1000 times
sills = []
ranges = []
grouped_variograms = {}
model = gs.Spherical()
for i in range(1000):
    random_points = extract_random_points(list(data), 1000)
    variogram = experimental_variogram(random_points)
    #print(variogram)
    plt.plot([distance for distance, gamma in variogram], [gamma for distance, gamma in variogram])
    plt.savefig("exp_variogram.png")
    plt.clf()
    for distance, gamma in variogram:
        if distance not in grouped_variograms:
            grouped_variograms[distance] = []
        grouped_variograms[distance].append(gamma)
    #calculating the mean variance for each distance
    mean_variances = []
    for distance, variances in grouped_variograms.items():
        mean_variances.append((distance, sum(variances) / len(variances)))

    #sorting the mean variances by distance
    mean_variances = sorted(mean_variances, key=lambda x: x[0])

    #plotting the new experimental variogram
    plt.plot([distance for distance, gamma in mean_variances], [gamma for distance, gamma in mean_variances])
    plt.savefig("mean_variogram.png")
    plt.clf()

    sill = max([gamma for _, gamma in variogram])
    range_ = 10000
    fit = fit_spherical_model(mean_variances, sill, range_)
    sills.append(sill)
    ranges.append(range_)
    #Saving the fit of the spherical model to a file
    with open("range_sill" + str(i) + ".dat", "w") as file:
        for distance, gamma in fit:
            file.write(str(distance) + " " + str(gamma) + "\n")
    #Saving the plot of the fit to an image
    plt.plot([distance for distance, _ in mean_variances], [gamma for _, gamma in mean_variances], label="Experimental Variogram")
    plt.plot([distance for distance, _ in fit], [gamma for _, gamma in fit], label="Spherical Model")
    plt.legend()
    plt.savefig("variogram" + str(i) + ".png")
    plt.clf()

#Saving the values of sill and range for each iteration
with open("sills.dat", "w") as file:
    for sill in sills:
        file.write(str(sill) + "\n")

with open("ranges.dat", "w") as file:
    for range_ in ranges:
        file.write(str(range_) + "\n")
