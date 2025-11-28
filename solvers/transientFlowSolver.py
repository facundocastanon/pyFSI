
from colorama import Fore, Back, Style

from pyFSI.solvers.solverBase import solverBase

# Solver for transient FSI simulations using the Newmark Integrator for the flow
class transientFlow(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)
        

    def solve(self):
        time = self._time
        
        # Integrate
        self.createFlowIntegrator()
        
        while time.value <= time.end:
            print("---> Solving the flow for time:", time.value)
            flowState = self.flowIntegrator.advance(time.span)
           
            print(Fore.BLUE +"---> Updating the flow")
            self._fsi.flow().updatePressure(flowState)
            
            self._fsi.flow().advance(time.span[1], flowState)
            
            time.advance()
            
            self._odb.write()
            
        self._odb.close()

