import numpy as np
import scipy.integrate as si

class scipyIVPSolid:
    def __init__(self, solid, time):
        self._solid = solid
        self._time = time
    
    def advance(self, tspan):

        fSol = si.solve_ivp(self._solid.RHS,
                            self._time.span,
                            self._solid.state,
                            method=self._solid.execution()['solver']['integrators']['solid']['type'],
                            atol=self._solid.execution()['solver']['integrators']['solid']['atol'],
                            rtol=self._solid.execution()['solver']['integrators']['solid']['rtol'])
        
        solidState = fSol.y[:, -1]
        
        # This line used to be be self._solid.rhs[self._solid.dof:self._solid.sof]
        # which is wrong, but it worked since there is a problem with the acceleration coupling
        # Therefore, I set the acceleration to zero, otherwise it may now work for large mass flows or large added mass (check the latter)
        acceleration = self._solid.rhs[self._solid.dof:self._solid.sof]
        
        return solidState, acceleration
    
    