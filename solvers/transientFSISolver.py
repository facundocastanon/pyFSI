import scipy.integrate as si
import numpy as np
import importlib
from pyFSI.execution.errors import *
from pyFSI.execution.utilities import subSection
from pyFSI.solvers.solverBase import solverBase
from colorama import Fore, Back, Style

# Solver for transient FSI simulations using the Newmark Integrator for the solid
class transientFSI(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)

    def solve(self):
        time = self._time
        ctrl = self.control
        
        self._fsi.setInitialConditions(time.span[1])
        
        # Choose the solver couping mode
        if ctrl['coupling']["type"] == 'implicit':
            self.itmax = ctrl['coupling'].get('maxIterations', 50)  # Default iterations
            print(Fore.BLUE + "---> Operating the transient FSI solver in implicit mode...")
        elif ctrl['coupling']["type"] == 'explicit':
            self.itmax = 1
            print(Fore.BLUE + "---> Operating the transient FSI solver in explicit mode...")
        else:
            errorExit("Unknown transient FSI solver coupling mode.")
            
        # Create the integrator objects
        self.createSolidIntegrator()
        self.createFlowIntegrator()
    
        # Integration loop
        while time.value <= time.end:
            subSection("Solving for time: " + str(time.value), where='high')

            self.iterate(time.span) 
            
            time.advance()
            
            self._odb.write()
            
        self._odb.close()
        

    def iterate(self, tspan):
        fsi = self._fsi
        #print(Fore.GREEN + "---> transientPyFSI solver iteration: " + str(iter))
        # 1) Solve the flow
        #print("---> Solving the flow")
        iter = 1
        lastIter = False
        while iter <= self.itmax:
            print(Fore.GREEN + "-> Starting iteration: " + str(iter))
            
            # Solve the flow
            print(Fore.BLUE + "---> Solving the flow...")
            flowState = self.flowIntegrator.advance(tspan)      
            formattedState = [f"{x:.3e}" for x in np.nditer(flowState)]
            print(Fore.BLUE + "--> Flow state is: ", formattedState)
            print(Fore.BLUE + "--> Updating the flow")
            flowConvergence = fsi.flow().updatePressure(flowState)
            
            # Solve the solid with the flow forcing (Block Gauss-Seidel)
            print(Fore.YELLOW + "---> Solving the solid...")
            solidState, solidAcc = self.solidIntegrator.advance(tspan)
            formattedState = [f"{x:.3e}" for x in np.nditer(solidState)]
            print(Fore.YELLOW + "--> Solid state is: ", formattedState)
            print(Fore.YELLOW + "--> Updating the solid")
            fsi.solid().update(tspan[1], solidState, solidAcc)

            # Update the regions with the solid kinematics
            fsi.flow().updateRegions()
            
            # Calculate the residuals at the interfaces  
            if self.control['coupling']['type'] == 'implicit':
                convergence = True
                for region in fsi.flow().regions:  
                    for boundary in region.boundaries():
                        if boundary.isCoupled():
                            for field in boundary.control['couple']:
                                if field == "Displacements":
                                    fnow = boundary.getDisplacements()
                                    fold = boundary.oldDisplacements()
                                elif field == "Velocities":
                                    fnow = boundary.getVelocities()
                                    fold = boundary.oldVelocities()                                
                                elif field == "Accelerations":
                                    fnow = boundary.getAccelerations()
                                    fold = boundary.oldAccelerations()                               
                                elif field == "Forces":
                                    fnow = boundary.getForces()
                                    fold = boundary.oldForces()                                  
                                
                                # Calculate the residual of the corresponding field
                                fres = self.residual(fnow, fold) 
                                
                                # Relax the fields
                                if lastIter == False:
                                    if 'relaxation' in self.control['coupling']:
                                            # print(Fore.BLUE + "---> Relaxing the fields...")
                                            if field in self.control['coupling']['relaxation']:
                                                boundary.relax(field, self.control['coupling']['relaxation'][field])
                                            
                                # Check convergence
                                if fres >= self.control['coupling']['tolerances'][field]:
                                    convergence = False
                                print(Fore.CYAN + "---> The residual for", field, "for region", region.name, "is:", format(fres, ".2e"))
            
            elif self.control['coupling']['type'] == 'explicit':
                lastIter = True
                                
            if lastIter == False:
                if iter == self.itmax or convergence == True:
                    if convergence == True:
                        subSection("The FSI system converged in " + str(iter) + " iterations." + "\n")
                    else:
                        subSection("The FSI system reached ITMAX. WARNING: Advancing anyway. " + "\n")
                    subSection("Performing the last iteration" + "\n")
                    lastIter = True
                    continue
                else:
                    iter += 1
                
            # If last iteration of converged, advance in time
            else:
                # Converged or ITMAX, advance in time
                fsi.flow().advance(tspan[1], flowState)
                fsi.solid().advance(solidState, solidAcc)
                break

                
            