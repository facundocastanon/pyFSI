# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#    p    #     version: 0.2
#    y    #     date: 1/11/2023
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   The class  works either as an analytical model or with a Calculix model
#   The model type is selected from the key "method" in the solid-solution dict
#   Supports reconstruction of displacements
# Warnings:
#   Modal damping is the same for every mode
#   Reference displacements fixed to L/10
#   Initial position of the beam is set to zero
# Optimize:
#   Delete info from dy and ddy
# --------------------------------------------------------------------------- #
import sys
import numpy as np
from numpy.core.numeric import ones_like
import scipy.optimize as opt
from scipy import interpolate
from scipy import optimize
import warnings
import scipy.integrate as si
import matplotlib.pyplot as plt
from colorama import Fore, Back, Style

from pyFSI.vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys
from pyFSI.models.solidModels.solidBase import solidModel
from pyFSI.post.prints import *

class spring(solidModel):
    def __repr__(self):
        return 'springModel'
    
    def __init__(self, execution, control, mesh, boundary, time):
        super().__init__(execution, control, mesh, boundary, time)  # Call the base class
        # ----- Public attributes ----- #
        
        self.type = 'spring'
        # Degrees of freedom
        self.dof = 1                 # Modal DOFs (one for the base motion)
        self.sof = 2 * self.dof      # State Space DOFS

        # Activate from zero if not defined
        if 'cutIn' not in self.control:
            self.control['cutIn'] = 0.0

        self.F = 0
        self.W = 0

        self.powers = range(len(self.control['coefficients']))
        
        self.varMap["forces"] = "F"
        self.varMap["power"] = "W"

    def force(self, displacement):
        # Calculate the spring force as a function of the displacement
        self.F = 0
        if self._time.value >= self.control['cutIn']: # Check the cutIn condition
            for i in self.powers:
                self.F += self.control['coefficients'][i] * displacement**(i+1)

        return self.F
