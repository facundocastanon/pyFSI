# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: eigen vector class
#    I    #     return: eigen vector object
# --------------------------------------------------------------------------- #
# Notes:
#   This class creates an eigenvector object inheriting from a numpy array
#   Note the this class inherits from field, which inherits from np.ndarray
#   Not strictly neccessary, but helps to understand inheritance (sorry)
# --------------------------------------------------------------------------- #
import numpy as np
import scipy.integrate as si

from pyFSI.vectors.field import field

class eigenVector(field):
    def __new__(cls, 
                nparray,
                calculate=False,
                mesh=None, 
                normalize=None, 
                mass=None, 
                factor=None, 
                first=None, 
                second=None, 
                third=None,
                fourth=None, 
                integralx=None, 
                integralL=None):

        # ----- Casting  ----- #
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        # Choose Normalization
        if normalize == 'length':
            normFactor = abs((mesh.L))**0.5
        elif normalize == 'mass':  # Typical mass normalization
            normFactor = (1 / si.simps(mass*nparray**2, mesh.x))**0.5
        elif normalize == 'general':  # Normalization for non-conventional BCs
            normFactor = factor
        else:
            normFactor = 1

        obj = np.asarray(nparray).view(cls) * normFactor

        # ----- Public attributes ----- #
        obj.normalization = normalize

        # If we choose to calculate, evaluate the derivatives. Note that we need to pass a mesh, otherwise we cannot evaluate  the derivatives and the integrals.
        if calculate:
            edgeOrder = 2
            # Choose between analytic or numerical derivatives and integrals
            # First derivative
            if first is not None:
                obj.d1 = normFactor * first
            else:
                obj.d1 = np.gradient(obj, mesh.x, edge_order=edgeOrder)
            # Second derivative
            if second is not None:
                obj.d2 = normFactor * second
            else:
                obj.d2 = np.gradient(obj.d1, mesh.x, edge_order=edgeOrder)
            # Third derivative
            if third is not None:
                obj.d3 = normFactor * third
            else:
                obj.d3 = np.gradient(obj.d2, mesh.x, edge_order=edgeOrder)
            # Fourth derivative
            if fourth is not None:
                obj.d4 = normFactor * fourth
            else:
                obj.d4 = np.gradient(obj.d3, mesh.x, edge_order=edgeOrder)
            # X Integral
            if integralx is not None:
                obj.ix = normFactor * integralx
            else:
                obj.ix = si.cumtrapz(obj, mesh.x, initial=0.0)
            # L Integral
            if integralL is not None:
                obj.iL = normFactor * integralL
            else:
                obj.iL = si.simps(obj, mesh.x)



        # ----- Private attributes ----- #
        obj._mesh = mesh

        return obj

    # Re-execute the integral after application of boundary conditions
    def correct(self):
        self.ix = si.cumtrapz(self, self._mesh.x, initial=0.0)
        self.iL = si.simps(self, self._mesh.x)

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return

    # Getters
    def mesh(self):
        return self._mesh
