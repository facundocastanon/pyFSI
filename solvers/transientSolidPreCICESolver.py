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
from pyFSI.math.ode import scipyIVPSolid

# Solver for transient FSI simulations
class transientSolidPreCICE(solverBase):
    def __init__(self, fsi, odb):
        super().__init__(fsi, odb)
        self.cleanCase()  # Clean the case for optimal linking

    def solve(self):
        interface = self._fsi.interface
        meshes = self._fsi.meshes
        time = self._time
        solid = self._fsi.participant
        
        # Create the solid integrator
        self.createSolidIntegrator()
        
        # preCICE defines timestep size of solver via precice-config.xml
        print(Fore.GREEN + "--> Initializing the preCICE interface for the solid ...")
        precice_dt = interface.initialize()
        
        
        if interface.is_action_required(prc.action_write_initial_data()):
            print(Fore.GREEN + "--> Writing available initial data  (solid side) ...")
            for mesh in meshes:
                mesh.writeFields(interface, solid)            
                interface.mark_action_fulfilled(prc.action_write_initial_data())

        interface.initialize_data()

        # Read initial data
        if interface.is_read_data_available():
            print(Fore.GREEN + "--> Reading available initial data  (solid side) ...")
            for mesh in meshes:
                mesh.readFields(interface)

        print(Fore.GREEN + "--> Setting the intial conditions (solid side) ...")
        totalTime = precice_dt  # Calculate the accumulated simluation time
        


        # # Write initial force data data
        # if interface.is_action_required(prc.action_write_initial_data()):
        #     print(Fore.GREEN + "--> Writing available initial data (solid side) ...")
        #     for mesh in meshes:
        #         mesh.writeFields(interface, solid)
        #     interface.mark_action_fulfilled(prc.action_write_initial_data())

        # Start the solution loop
        print(Fore.GREEN + "--> Starting the preCICE loop (solid side) ...")
        while interface.is_coupling_ongoing():
            # When an implicit coupling scheme is used, checkpointing is required
            if interface.is_action_required(prc.action_write_iteration_checkpoint()):
                interface.mark_action_fulfilled(prc.action_write_iteration_checkpoint())

            # Solve the flow
            print(Fore.GREEN + "--> Calling the solver (solid side) ...")
            solidState, accState = self.solidIntegrator.advance(precice_dt) 
            solid.update(totalTime, solidState, accState)     
        
            # Write the force in the precice interface
            print(Fore.GREEN + "--> Writing the solid kinematics to the partner mesh (solid side) ...")
            for mesh in meshes:
                mesh.writeFields(interface, solid)

            # Advance in time
            print(Fore.GREEN + "--> Advancing the precice interface (solid side) ...")
            precice_dt = interface.advance(precice_dt)  # Execute the mapping

            # Read the displacement and velocities and write the the pyFSI boundary objects
            print(Fore.GREEN + "--> Reading the solid fields (solid side) ...")
            for mesh in meshes:
                mesh.readFields(interface)


            # Check convergence and confirm advancing
            print(Fore.GREEN + "--> Cheking convergence (solid side) ...")
            if interface.is_action_required(prc.action_read_iteration_checkpoint()):  # i.e. not yet converged
                print(Fore.GREEN + '--> The implicit preCICE step has NOT converged (solid side)')
                interface.mark_action_fulfilled(prc.action_read_iteration_checkpoint())
            else:  # converged, timestep complete
                print(Fore.GREEN + '--> The implicit preCICE step has converged (solid side)')
                totalTime += precice_dt
                solid.advance(solidState, accState)                
                time.advance(dt=precice_dt)
                self._odb.write()

        print(Fore.GREEN + "Exiting precice...")

        interface.finalize()

        self._odb.close()

        return solid

    # Clean the case (neccessary for optimal linking betwen the participants)
    def cleanCase(self):
        print(Fore.GREEN + "--> Allcleaning preCICE...")
        for p in pathlib.Path(".").glob("P*.log"):
            p.unlink()
        try:
            shutil.rmtree("precice-run")
        except OSError as e:
            print(Fore.GREEN + "Error: %s - %s." % (e.filename, e.strerror))
