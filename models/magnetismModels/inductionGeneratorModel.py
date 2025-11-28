import numpy as np
from pyFSI.post.prints import printArray
from pyFSI.models.magnetismModels.magnetismBase import magnetismModel
from pyFSI.util.misc import *

class inductionGenerator(magnetismModel):
    """Electromagnetic induction generator class."""
    def __init__(self, execution, control, mesh, time):
        """Creates and Electromagnetic Induction object. 

        Args:
            execution (dictionary): [description]
            control (dictionary): [description]
            mesh (fsiMesh1D): [description]
            time (time): [description]
        """
        super().__init__(execution, control, mesh, time)

        # ----- Public attributes ----- #
        self.V = 0.0  # Voltage
        self.W = 0.0  # Power
        
        # ----- Private attributes ----- #       
        C = self.control

        # Induction coefficient and resistance
        self._indu = C['dphi'] * C['area'] * C['loops']
        self._resi = self.control['resistance']
        
        
        # Scale the inductance and the resistance
        if "scale" in C:
            self._indu *= C['scale']
            self._resi *= C['scale']

        # Output variable mapping
        self.varMap["voltage"] = "V"
        self.varMap["power"] = "W"

        self._coupled = getBool(C, "coupled")

    def update(self, vel):
        """Updates the power and voltage
        Args:
            vel (ndarray): Velocity of the associated solid object
        """
        # Test scaling the resistance   
        # self._resi = 100 + 3 * self._time.value
        # print(self._resi)
        self.V =  self._indu * vel
        self.W = self.V**2 / self._resi


    def beamConcentratedForce(self, pos, vel):
        """Creates a beam force vector that acts as a point force.

        Args:
            pos (ndarray): Current position of the beam
            vel (ndarray): Current velocity of the beam

        Returns:
            ndarray: Force vector
        """

        # Update the power and the voltage
        self.update(vel)

        # Calculate the mechanical viscous like force (avoid devide by vel)
        pointForce = self._indu**2 * vel / self._resi 
        
        # Create the force vector, add the force. Expensive but simple.
        if self._coupled:
            force = -pointForce  # Negative, it opposes the motion.
        else:
            force = 0.0
        
        return force

