
from colorama import Fore, Back, Style
import numpy as np
from pyFSI.solvers.solverBase import solverBase
from pyFSI.math.ode import generalizedAlpha

# Solver for transient FSI simulations using the Newmark Integrator for the solid
class transientSolid(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)
        

    def solve(self):
        time = self._time
        
        # Integrate
        self.createSolidIntegrator()
                
        while time.value <= time.end:
            print("---> Solving the structure for time:", time.value)
            solidState, solidAcc = self.solidIntegrator.advance(time.span)
           
            formattedState = [f"{x:.3e}" for x in np.nditer(solidState)]
            print(Fore.YELLOW + "--> Solid state is: ", formattedState)
            
            self._fsi.solid().update(time.span[1], solidState, solidAcc)
            
            self._fsi.solid().advance(solidState, solidAcc)
            
            time.advance()
            
            self._odb.write()
            
        self._odb.close()

