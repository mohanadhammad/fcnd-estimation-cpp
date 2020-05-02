import numpy as np


def parseData(filename):
    file = open(filename,'r')
    lines = file.readlines()

    samples = []
    isFirstLine = True

    for line in lines:
        if isFirstLine is False:
            line = line.strip()
            args = line.split(',')
            samples.append(float(args[1]))
        isFirstLine = False

    return samples

def calcMeanAndVariance(samples):
    mu = 0.
    for i in range(len(samples)):
        mu += samples[i]
    mu /= len(samples)

    var = 0.
    for i in range(len(samples)):
        var = var + (samples[i] - mu)**2
    var /= len(samples)

    return mu, var

if __name__ == "__main__":
    gpsData = parseData('../config/log/Graph1.txt')
    imuData = parseData('../config/log/Graph2.txt')

    gpsMean, gpsVar = calcMeanAndVariance(gpsData)
    imuMean, imuVar = calcMeanAndVariance(imuData)

    print("np.sqrt(gpsVar) = ", np.sqrt(gpsVar))
    print("np.sqrt(imuVar) = ", np.sqrt(imuVar))

