import numpy as np
import scipy.integrate as si

class scipyIVP:
    def __init__(self, solid, time):
        self._solid = solid
        self._time = time
    
    def advance(self, tspan):

        fSol = si.solve_ivp(self._solid.RHS,
                            tspan,
                            self._solid.state,
                            method=self._solid.execution()['solver']['integrator'][1],
                            atol=self._solid.execution()['solver']['atol'],
                            rtol=self._solid.execution()['solver']['rtol'])
        
        solidState = fSol.y[:, -1]
            
        return solidState, self._solid.rhs[self._solid.dof:self._solid.sof]