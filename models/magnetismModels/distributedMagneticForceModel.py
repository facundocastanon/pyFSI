
import numpy as np
from pyFSI.models.magnetismModels.magnetismBase import magnetismModel

class distributedMagneticForce(magnetismModel):

    def __init__(self, execution, control, mesh, time):
        super().__init__(execution, control, mesh, time)  # Call the base class

        # ----- Private attributes ----- #
        self.k1 = control['k1']
        self.k2 = control['k2']
        self.k3 = control['k3']

    def beamDistributedForce(self, beam):
        y = beam.y['mid']
        force = self.k1 * y + self.k2 * y**2 + self.k3 * y**3
        return force

