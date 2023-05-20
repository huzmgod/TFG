import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import random
import math
import gstools as gs
from collections import Counter

# Other required distributions
distributions = [stats.norm, stats.gamma, stats.beta, stats.expon, stats.chi2, stats.weibull_min]

# Load data
data = np.loadtxt("datafiles/vel_interpolada_DEM_ordenada.dat")

def euclideanDistance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

def extractRandomPoints(data, n):
    return random.sample(list(data), n)

def experimentalVariogram(points):
    variogram = []
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            x1, y1, v1 = points[i]
            x2, y2, v2 = points[j]
            distance = euclideanDistance(x1, y1, x2, y2)
            if 50 <= distance <= 50000:
                variogram.append((distance, 0.5 * (v1 - v2)**2))
    return sorted(variogram, key=lambda x: x[0])

def sphericalModel(distance, sill, range_, nugget=0):
    if distance <= range_:
        return nugget + (sill - nugget) * ((3 * distance)/(2 * range_) - (distance / range_)**3)
    else:
        return sill

def fitSphericalModel(variogram, sill, range_, nugget=0):
    fit = []
    for distance, gamma in variogram:
        fit.append((distance, sphericalModel(distance, sill, range_, nugget)))
    return fit

sills = []
ranges = []

for i in range(1000):
    print(f"Iteration {i + 1} of 1000")
    randomPoints = extractRandomPoints(data, 1000)
    variogram = experimentalVariogram(randomPoints)
    sill = max(gamma for _, gamma in variogram)
    range_ = 10000
    fit = fitSphericalModel(variogram, sill, range_)
    sills.append(sill)
    ranges.append(range_)

def plotHistogramAndFit(data, bins, distribution, params, title):
    plt.hist(data, bins, alpha=0.6, color='g', density=True)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    fittedData = distribution.pdf(x, *params[:-2], loc=params[-2], scale=params[-1])
    plt.plot(x, fittedData, 'r-', lw=2)
    plt.title(title)
    plt.savefig(f"{title.replace(' ', '_')}.png")
    plt.clf()

def fitDistribution(data, distribution):
    params = distribution.fit(data)
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Calculate the PDF (fitted distribution)
    pdfFitted = distribution.pdf(np.linspace(min(data), max(data)), *arg, loc=loc, scale=scale)

    return params, np.sum((pdfFitted - data)**2)

def findBestDistribution(data):
    bestDistribution = None
    bestParams = None
    smallestError = float('inf')

    for distribution in distributions:
        params, error = fitDistribution(data, distribution)
        if error < smallestError:
            smallestError = error
            bestDistribution = distribution
            bestParams = params

    return bestDistribution, bestParams

bins = range(10, 101)

# Find best fit distribution for sills
print("Finding best fit for sills...")
bestDistributionSills, bestParamsSills = findBestDistribution(sills)
plotHistogramAndFit(sills, bins, bestDistributionSills, bestParamsSills, "Sills Histogram and Best Fit")

# Find best fit distribution for ranges
print("Finding best fit for ranges...")
bestDistributionRanges, bestParamsRanges = findBestDistribution(ranges)
plotHistogramAndFit(ranges, bins, bestDistributionRanges, bestParamsRanges, "Ranges Histogram and Best Fit")

print(f"Sills best fit: {bestDistributionSills.name}, parameters: {bestParamsSills}")
print(f"Ranges best fit: {bestDistributionRanges.name}, parameters: {bestParamsRanges}")
