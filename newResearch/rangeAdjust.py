import random
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import euclidean
from scipy.optimize import minimize
from scipy.stats import norm, gamma, beta
import gstools as gs
from sklearn.mixture import GaussianMixture

# Load data
data = np.loadtxt("datafiles/vel_interpolada_DEM_ordenada.dat")

# Calculating the euclidean distance between the points
def euclidean_distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

# Extracting 1000 random points from the data
def extract_random_points(data, n):
    return random.sample(data.tolist(), n)

# Calculating the experimental variogram
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

# Fitting the experimental variogram to a spherical theoretical model
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

# Repeating the process 1000 times
sills = []
ranges = []
for i in range(1000):
    random_points = extract_random_points(data, 1000)
    variogram = experimental_variogram(random_points)

    sill = max([gamma for _, gamma in variogram])
    range_ = 10000
    fit = fit_spherical_model(variogram, sill, range_)

    sills.append(sill)
    ranges.append(range_)

# After obtaining all the sill and range values, plot a histogram for each
# This histogram will represent the distribution of the data
sns.histplot(sills, bins=50, kde=True)
plt.xlabel('Sill values')
plt.ylabel('Frequency')
plt.title('Histogram of Sill Values')
plt.savefig("sill_histogram.png")
plt.clf()

sns.histplot(ranges, bins=50, kde=True)
plt.xlabel('Range values')
plt.ylabel('Frequency')
plt.title('Histogram of Range Values')
plt.savefig("range_histogram.png")
plt.clf()

# Next, we will try to fit a known distribution to these data
# We will try Normal, Gamma, and Beta distributions and choose the one with the least error
distributions = [norm, gamma, beta]

def fit_distribution(data, dist):
    params = dist.fit(data)
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]
    
    # Calculate the PDF (fitted distribution)
    pdf_fitted = dist.pdf(np.linspace(min(data), max(data)), *arg, loc=loc, scale=scale)
    
    return params, np.sum((pdf_fitted - data)**2)

def find_best_distribution(data):
    best_distribution = None
    best_params = None
    best_error = np.inf
    
    for dist in distributions:
        params, error = fit_distribution(data, dist)
        if error < best_error:
            best_distribution = dist
            best_params = params
            best_error = error
            
    return best_distribution, best_params

sill_distribution, sill_params = find_best_distribution(sills)
range_distribution, range_params = find_best_distribution(ranges)

print(f"Sill Distribution: {sill_distribution.name}, Parameters: {sill_params}")
print(f"Range Distribution: {range_distribution.name}, Parameters: {range_params}")
