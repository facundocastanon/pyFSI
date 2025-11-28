
import numpy as np
from pyFSI.models.magnetismModels.magnetismBase import magnetismModel
from pyFSI.post.prints import *


class distributedMagneticStiffness(magnetismModel):

    def __init__(self, execution, control, mesh, time):
        super().__init__(execution, control, mesh, time)  # Call the base class

        # ----- Private attributes ----- #
        self.k1 = control['k1']

    def beamDistributedForce(self, beam):
        # This complex logic to modify the eigenvalue avoids the use of beta
        EI = beam.control["EI"]
        m = beam.control["m"]
        for i, value in enumerate(beam.eigen.values):
            w2 = value**2  # Original squared omega
            wk2 = (self.k1 / EI) * ((EI) / m) # Additional omega squared
            w2 += wk2 # New omega squared
            if w2 > 0:
                dwk = w2**0.5 - value   # Delta Omega
                value += dwk # New Omega
                np.put(beam.eigen.values,[i],[value])
            else:
                w = np.emath.sqrt(w2)
                np.put(beam.eigen.values,[i],[w])
        
        # We return 0 because there is no net force but a modification of the
        # eigenvalue of the system
        return 0

