# --------------------------------------------------------------------------- #
#    p    #     version: 0.3.0
#    y    #     date: 07/12/2020
#    F    #     author: Martin Saravia
#    S    #     description: preCICE solver
#    I    #     return: solver object
# --------------------------------------------------------------------------- #
# Notes:
#   Staggered solver for fsi using preCICE and Calculix
#   a) The solid is solved using Calculix
# --------------------------------------------------------------------------- #

import importlib, shutil, os, pathlib
import scipy.integrate as si
import numpy as np
import precice as prc
from colorama import Fore, Back, Style

from pyFSI.solvers.solverBase import solverBase


# Solver for transient FSI simulations
class transientFlowPreCICE(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)
        self.cleanCase()  # Clean the case for optimal linking

    def solve(self):
        iteration = 0
        interface = self._fsi.interface
        meshes = self._fsi.meshes
        time = self._time
        flow = self._fsi.participant
        
        # Create the integrator
        self.createFlowIntegrator()
        
        # preCICE defines timestep size of solver via precice-config.xml
        print(Fore.GREEN + "--> Initializing the preCICE interface for the flow...")
        precice_dt = interface.initialize()
        
        # Write initial force data data
        if interface.is_action_required(prc.action_write_initial_data()):
            print(Fore.GREEN + "--> Writing available initial data ...")
            for mesh in meshes:
                mesh.writeFields(interface, flow)
            interface.mark_action_fulfilled(prc.action_write_initial_data())
            
        interface.initialize_data()  
        
        # Read initial data
        if interface.is_read_data_available():
            print(Fore.GREEN + "Reading available initial data  (fluid side)...")
            for mesh in meshes:
                mesh.readFields(interface)
  

        print(Fore.GREEN + "--> Setting the intial conditions (fluid side)...")
        totalTime = precice_dt  # Total simluation time


        # Start the solution loop
        print(Fore.GREEN + "--> Starting the preCICE loop (fluid side)...")
        while interface.is_coupling_ongoing():
            # When an implicit coupling scheme is used, checkpointing is required
            if interface.is_action_required(prc.action_write_iteration_checkpoint()):
                interface.mark_action_fulfilled(prc.action_write_iteration_checkpoint())
                
            # Solve the flow
            flow.updateRegions()
            flowState = self.flowIntegrator.advance(self._time.span)

            # Update the flow state
            print(Fore.GREEN + "--> Updating the flow object (fluid side)...")
            flowConvergence = flow.updatePressure(flowState)

            # Write the force in the precice interface
            print(Fore.GREEN + "--> Writing the fluid force to the partner mesh (fluid side)...")
            for mesh in meshes:
                mesh.writeFields(interface, flow)

            # Advance in time
            print(Fore.GREEN + "--> Advancing the precice interface (fluid side)...")
            precice_dt = interface.advance(precice_dt)  # Execute the mapping

            # Read the displacement and velocities and write the the pyFSI boundary objects
            print(Fore.GREEN + "--> Reading the solid fields (fluid side)...")
            for mesh in meshes:
                mesh.readFields(interface)

            # # Update the region geometry after updating the boundaries
            # print(Fore.GREEN + "--> Updating the region dynamics...")
            # for i, region in enumerate(self.fsi.flow().regions):
            #     region.update()

            # Check convergence and confirm advancing
            print(Fore.GREEN + "--> Cheking convergence (fluid side)...")
            if interface.is_action_required(prc.action_read_iteration_checkpoint()):  # i.e. not yet converged
                interface.mark_action_fulfilled(prc.action_read_iteration_checkpoint())
                print(Fore.GREEN + '--> The implicit preCICE step has NOT converged (flow side). Iteration: ' + str(iteration))
                iteration += 1
            else:  # converged, timestep complete
                print(Fore.GREEN + '--> The implicit preCICE step has converged (flow side)')
                totalTime += precice_dt
                flow.advance(self._time.span[1], flowState)              
                time.advance(dt=precice_dt)
                self._odb.write()
                iteration = 0

        print(Fore.GREEN + "Exiting precice...")

        interface.finalize()

        self._odb.close()

        return flow

    # Clean the case (neccessary for optimal linking betwen the participants)
    def cleanCase(self):
        print(Fore.GREEN + "--> Allcleaning preCICE...")
        for p in pathlib.Path(".").glob("P*.log"):
            p.unlink()
        try:
            shutil.rmtree("precice-run")
        except OSError as e:
            print(Fore.GREEN + "Error: %s - %s." % (e.filename, e.strerror))
