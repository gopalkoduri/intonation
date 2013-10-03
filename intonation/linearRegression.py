from __future__ import division
import numpy as np
from numpy.linalg.linalg import pinv

def computeCost(X, y, theta):
    """
    X: matrix of features, one sample per row (with bias unit)
    """
    numSamples = len(y)
    error = X*theta - y
    cost = sum(np.multiply(error, error))/(2*numSamples)
    return cost


def gradientDescent(X, y, theta, alpha, maxIter=100):
    """
    X: matrix of features, one sample per row (without bias unit)
    y: values (continuous) corresponding to rows (samples) in X
    theta: parameters as row vector/matrix
    alpha: learning rate
    iterations: num of iterations

    returns [theta, costHistory]
    """
    print "Note: Remember that features should be normalized for gradient\
    descent to work well..."
    numSamples = y.size
    X = np.insert(X, 0, np.ones(numSamples), axis=1)

    for i in xrange(maxIter):
        newTheta = theta - (alpha/numSamples)*X.T*(X*theta - y)
        if (theta == newTheta).all():
            print i
            break

        theta = newTheta

    cost = computeCost(X, y, theta)
    return [theta, cost]

def normalEquation(X, y):
    """
    X: matrix of features, one sample per row (without bias unit)
    y: values (continuous) corresponding to rows (samples) in X
    """
    numSamples = y.size
    X = np.insert(X, 0, np.ones(numSamples), axis=1)

    return pinv(X)*y

def predict(X, theta):
    return X*theta

