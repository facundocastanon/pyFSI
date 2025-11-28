import numpy as np
from pyFSI.vectors.eigen import eigenSystemVector as es
from pyFSI.solvers.solverBase import solverBase
from pyFSI.post.prints import *

class eigen(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)

    def solve(self):
        time = self._time
        fsi = self._fsi
        oreg = fsi.solid().mesh().OREG
        while time.value <= time.end:
            for obj in oreg.objects:
                if 'parameter' in obj.control:
                    obj.parametrize(time)
            print("  Solving the MFSI case:", fsi.name, "for time ", time.value)
            fsi.update()
            fsi.solve(self)
            time.advance()
            self._odb.write()

        self._odb.close()

    def response(self, x0, time, modes=None):
        """ This function returns the evolution of the FSI system as a sum 
        of harmonic functions for a unit initial condition in all modes.

        Args:
            x0 (array): Initial condition
            time (array): Time span to calculate
            modes (list): Modes to use in the response calculation

        Returns:
            Array: Evolution of the FSI modal response
        """
        fsi = self._fsi
        
        # Select which modes are included in the responsa calculationn
        if modes:
            print("---> Linear FSI response will be calculated with modes: ", modes)
            pass
        else:
            modes = range(fsi.ES.size)
        
        # Initialize the transient generalized coordinates (all FSI modes)   
        xt = np.zeros((fsi.ES.size, len(time)), dtype='complex')
        
        # Calculate the generalized coordinate for the selected modes
        ai = np.zeros(fsi.ES.size)
        for i in modes:
            ai[i] = np.dot(fsi.Y.T[i, :], x0)
            
        # Calculate the response using the selected modes
        for t, tval in enumerate(time):
            for i in modes:
                xt[:, t] += fsi.X[:, i] * np.exp(fsi.R[i] * tval) * ai[i]
                print("ai: ", ai)
                print("mode: ", i)
                print("eigenvalue", fsi.R[i])
                print("Eigenvector", fsi.X[:,i])
                
        return xt
    
  