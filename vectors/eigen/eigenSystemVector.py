import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
from pyFSI.vectors.field import field
from pyFSI.vectors.eigen.eigenValue import *
from pyFSI.vectors.eigen.eigenVector import *

from abc import ABCMeta, abstractmethod


class eigenSystemVector(object):

    def __repr__(self):
        return 'eigenSystemVector'

    def __init__(self, evalues, evectors, name=None, calculate=False,           complex=False):

        # ----- Public attributes ----- #
        self.name = name
        self.size = len(evalues)
        self.values = np.array(evalues)  # Eigenvalues
        self.vectors = np.zeros(self.size, dtype=object) # Eigenvectors
        for i in range(self.size):
            self.vectors[i] = evectors[i]
        # The system has complex eigenvalues and eigenvectors
        if complex:
            self.valuesReal = np.zeros(self.size)
            self.valuesImag = np.zeros(self.size)
            self.vectorsReal = np.zeros(self.size, dtype=object) # Eigenvectors
            self.vectorsImag = np.zeros(self.size, dtype=object) # Eigenvectors
            for i in range(self.size):
                self.valuesReal[i] = self.values[i].real
                self.valuesImag[i] = self.values[i].imag
                self.vectorsReal[i] = self.vectors[i].real
                self.vectorsImag[i] = self.vectors[i].imag
            
            
        if calculate:
            self.d1 = np.zeros(self.size, dtype=object)  # First derivative
            self.d2 = np.zeros(self.size, dtype=object)
            self.d3 = np.zeros(self.size, dtype=object)
            self.d4 = np.zeros(self.size, dtype=object)
            self.ix = np.zeros(self.size, dtype=object)  # indefinvaite integrals
            self.iL = np.zeros(self.size)  # definite integrals in L
            self.N = np.zeros((self.size, self.size), dtype=object)  # Norm
            # Assemble the eigensystem and its derivatives
            for i in range(self.size):
                self.vectors[i] = evectors[i]
                self.d1[i] = evectors[i].d1
                self.d2[i] = evectors[i].d2
                self.d3[i] = evectors[i].d3
                self.d4[i] = evectors[i].d4
                self.ix[i] = evectors[i].ix
                self.iL[i] = evectors[i].iL
        
        # ----- Private attributes ----- #
        self._mesh = evectors[0].mesh()
        
    # This method allows us to create the eigensystem from the solution given by numpy.eig. Note that we still do not allow to normalize.
    @classmethod
    def fromNumpy(cls, wArray, vArray):
        evaluesList = []   # Create the list of eivgenvalue objects
        evectorsList = []
        complex = False
        for i, value in enumerate(wArray):
            evaluesList.append(eigenValue(value)) # Append the object
            # Eigenvectors are the columns of the vArray
            evectorsList.append(eigenVector(vArray[:, i])) 
            if np.iscomplex(value):
                complex = True
        return cls(evaluesList, evectorsList, complex=complex)                    


    def correct(self):  # Correct the integral after setting the BCs
        for i in range(self.size):
            self.vectors[i].correct()
            self.N[i, i] = si.simps(self.vectors[i]**2, self._mesh.x)  # OrthoNorm

    # Reconstruct displacements from the state vector and the modes
    def reconstruct(self, state, modes=None):
        # Get the size of the mesh (number of elements in the eigenvectors)
        size = self.vectors[0].size 
        
        # Create a field of zeros
        fx = field(np.zeros(size))  
        
        # Select modes used in the reconstruction
        if modes is None:
            modes = range(self.size)
            
        # Reconstruct the deformed shape using mode superposition
        # This is slow, but easy to read (vectorize for speed)
        for i in modes:
            fx += state[i] * self.vectors[i]
        
        return fx

    def reconstructDdx(self, state):
        size = self.vectors[0].size  # Number of elements of eigenvectors
        dx = field(np.zeros(size))
        for i in range(self.size):
            dx += state[i] * self.d2[i]
        self.dx = dx
        return dx
    
    def plot(self):
        fig, ax = plt.subplots()
        legend = []
        for i, vector in enumerate(self.vectors):
            legend.append("Mode " + str(i)  + ": " + 
                          str(self.values[i].hz().round(2)) + " Hz")
            ax.plot(self._mesh.x, vector)
        ax.legend(legend)
        ax.set_title("Eigenvectors")
        plt.show()

    # Setters for boundary conditions
    def setV(self, rg, mode, val):
        self.vectors[mode][rg] = val

    def setD1(self, rg, mode, val):
        self.d1[mode][rg] = val

    def setD2(self, rg, mode, val):
        self.d2[mode][rg] = val

    def setD3(self, rg, mode, val):
        self.d3[mode][rg] = val

    def setD4(self, rg, mode, val):
        self.d4[mode][rg] = val



