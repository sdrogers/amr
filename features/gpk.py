"""Kernel functions that calculate the covariance between two inputs."""

import numpy as np
import abc

from scipy.spatial import distance


class BaseKernel(abc.ABC):

    """ A Gaussian Process kernel.
       Attributes:
           hypers (list)
    """

    @abc.abstractmethod
    def __init__(self):
        """ Create a GPKernel. """
        self._n_hypers = 0

    @abc.abstractmethod
    def cov(self, X1, X2, hypers=None):
        """ Calculate the covariance. """
        return np.zeros((len(X1), len(X2)))

    @abc.abstractmethod
    def fit(self, X):
        return self._n_hypers

