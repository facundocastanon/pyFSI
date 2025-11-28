# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#    p    #     version: 0.4
#    y    #     date: 20/11/2023
#    F    #     author: Martin Saravia
#    S    #     description: Damper Blass
#    I    #     return: damperModel object
# --------------------------------------------------------------------------- #
# Notes:
#   This class represents a polynomial damper with general coefficients
#
# Warnings:
#   
# Optimize:
#   It should be updated automatically by the solver and connected via boundaries
#   It is not good to update the power when other object calls self.force.
# --------------------------------------------------------------------------- #
import numpy as np
from colorama import Fore, Back, Style
from pyFSI.models.solidModels.solidBase import solidModel
from pyFSI.post.prints import *

class damper(solidModel):
    def __repr__(self):
        return 'damperModel'
    
    def __init__(self, execution, control, mesh, boundary, time):
        super().__init__(execution, control, mesh, boundary, time)  # Call the base class
        # ----- Public attributes ----- #

        self.type = 'damper'
        
        # Degrees of freedom
        self.dof = 1                 # Modal DOFs (one for the base motion)
        self.sof = 2 * self.dof      # State Space DOFS

        # Force and power in the damper
        self.F = 0
        self.W = 0

        # Activate from zero if not defined
        if 'cutIn' not in self.control:
            self.control['cutIn'] = 0.0

        self.powers = range(len(self.control['coefficients']) )
        print(self.control['coefficients'])
        
        self.varMap["forces"] = "F"
        self.varMap["power"] = "W"


    def force(self, velocity):
        # Calculate the damper force as a function of the displacement
        self.F = 0
        if self._time.value >= self.control['cutIn']:
            for i in self.powers:
                self.F += self.control['coefficients'][i] * velocity**(i+1)

        # Calculate the dissipated power
        self.W = self.F * velocity

        print(" power in the damper: ", self.W)

        return self.F
