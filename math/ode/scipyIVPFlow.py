import numpy as np
import scipy.integrate as si

class scipyIVPFlow:
    def __init__(self, flow, time):
        self._flow = flow
        self._time = time
    
    def advance(self, tspan):

        # Initialize flow parameters that are constant through the iteration
        self._flow.initialize()
        
        # Integrate the ODE
        fSol = si.solve_ivp(self._flow.RHS,
                            self._time.span,
                            self._flow.state,
                            method=self._flow.execution()['solver']['integrators']['flow']['type'],
                            atol=self._flow.execution()['solver']['integrators']['flow']['atol'],
                            rtol=self._flow.execution()['solver']['integrators']['flow']['rtol'])
        
        flowState = fSol.y[:, -1]
        
        return flowState
    
    