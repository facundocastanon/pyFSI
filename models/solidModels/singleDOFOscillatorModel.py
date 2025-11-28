# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#    p    #     version: 0.2
#    y    #     date: 03/11/2021
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   The class  works either as an analytical model or with a Calculix model
#   The model type is selected from the key "method" in the solid-solution dict
#   Supports reconstruction of displacements
# Warnings:
#   Modal damping is the same for every mode
#   Reference displacements fixed to L/10
#   Initial position of the beam is set to zero
# Optimize:
#   Delete info from dy and ddy
# --------------------------------------------------------------------------- #
import sys
import numpy as np
from numpy.core.numeric import ones_like
import scipy.optimize as opt
from scipy import interpolate
from scipy import optimize
import warnings
import scipy.integrate as si
import matplotlib.pyplot as plt
from colorama import Fore, Back, Style

from pyFSI.vectors.eigen import eigenValue as eval, eigenVector as evec, eigenSystemVector as esys
from pyFSI.models.solidModels.solidBase import solidModel
from pyFSI.post.prints import *

class singleDOFOscillator(solidModel):
    def __repr__(self):
        return 'singleDOFOscillatorModel'
    
    def __init__(self, execution, control, mesh, boundary, time):
        super().__init__(execution, control, mesh, boundary, time)  # Call the base class
        # ----- Public attributes ----- #
        # Degrees of freedom
        self.dof = 1     # Modal DOFs (one for the base motion)
        self.sof = 2 * self.dof                # State Space DOFS
        self.eigen = None                      # Eigensystem
        self.betas = None

        # Modal Functions
        self.k = np.zeros(self.dof, dtype=object)   # Stiffness vector
        self.c = np.zeros(self.dof, dtype=object)   # Damping vector
        self.m = np.zeros(self.dof, dtype=object)   # Mass vector

        # Galerkin System
        self.K = np.zeros((self.dof, self.dof))     # Modal Mass Matrix
        self.C = np.zeros((self.dof, self.dof))     # Modal Stiffness Matrix
        self.M = np.zeros((self.dof, self.dof))     # Modal Damping Matrix
        self.F = np.zeros(self.dof)                 # Modal Force
        
        # Modal State Space System
        self.S = np.zeros((self.sof, self.sof))     # State Space Matrix

        # State Space variables
        self.a = None                       # Generalized displacement
        self.da = None                      # Generalized velocity
        self.dda = None                     # Generalized acceleration

        # Physical variables
        # self.xi = None      # Initial position
        # self.xf = None      # Final position
        self.y = {}         # Current position of top, bot and mid surfaces
        self.y0 = {}        # Initial position
        self.dy = {}        # Current velocity of top, bot and mid surfaces
        self.ddy = {}       # Current accelerations of top, bot and mid surfaces
        self.v = {}         # Displacement vector         

        # Energy and power (elastic, kinetic, gravity)
        self.energy = None
        self.power = None
        
        # Dimensional analysis
        self.lRef = None  # Reference length
        self.tRef = None  # Reference thickness
        self.uRef = None  # Reference something
        self.vRef = None  # Reference velocity

        # Output variable mapping
        self.varMap["mesh"] = "_mesh.x"
        self.varMap["eigenValues"] = "eigen.values"
        self.varMap["eigenVectors"] =  "eigen.vectors"
        self.varMap["displacements"] = "y['mid']"
        self.varMap["velocities"] = "dy['mid']"
        self.varMap["accelerations"] = "ddy['mid']"
        self.varMap["numbers"] = "dimNumbers"
        self.varMap["energy"] = "energy"
        self.varMap["power"] = "power"
        self.varMap["state"] = "state"
        
        # These variables are written only at the beggining of the analysis 
        self.varMap["initial"] = ["eigenVectors", "eigenValues", "mesh"]

        # ----- Private attributes ----- #
        # Flags
        self._updated = False
        # Messages
        if self._debug:
            print("     WARNING: System initial conditions set to zero. ")

        # ----- Procedures ----- #
        # Add calculated parameters to the control dictionary


        # Create the eigensystem (it doet not use boundary info, just the mesh)
        self._createEigenSystem(control)
        
        # Calculate the dimensionless numbers
        self.calcNumbers()
        
        # Calculate the state matrix S
        self._createStateSpaceMatrix()  
        
        # Calculate vectors k, m, c and F for linear stability analysis
        self._createModalFunctions()  
        
        # Set the initial conditions and the boundary
        self.setInitialConditions()

    # ----- PUBLIC METHODS ---

    # Set the initial conditions
    def setInitialConditions(self, t0=0):
        # Initial state conditions
        initialState = np.zeros(self.sof)
        self.state = initialState
        self.a = initialState[0:self.dof]
        self.da = initialState[self.dof:self.sof]
        self.dda = np.zeros(self.dof)
        
        # Initial energy state
        self.energy0 = np.zeros(4)
        
        # Calculate the initial position
        self.y0['mid'] = self.eigen.reconstruct(initialState)
        
        # Update the solid
        self.update(None, initialState, self.dda)
        

    # Update the kinematics without updating the state (no convergence)
    def update(self, time, state, acc):
        x = self._mesh.x
        C = self.control
 
        # Update the state variables
        a = state[0:self.dof]
        da = state[self.dof:self.sof]
        dda = acc
        
        # Update the positions, displacements, velocities and accelerations
        self.y['mid'] = self.eigen.reconstruct(a)
        self.v['mid'] = self.y['mid'] - self.y0['mid']
        self.dy['mid'] = self.eigen.reconstruct(da)
        self.ddy['mid'] = self.eigen.reconstruct(dda)
 
        # Update the second spatial derivative (avoid reconstruction)
        ddx = self.eigen.reconstructDdx(a)
        
        # Update the associated boundaries
        for boundary in self.boundary:
            boundary.update("Displacements", self.v['mid'])
            boundary.update("Velocities", self.dy['mid'])
            boundary.update("Accelerations", self.ddy['mid'])
             
    def advance(self, state, acc):
        # Converged, so update the state variables
        self.state = state
        self.a = state[0:self.dof]
        self.da = state[self.dof:self.sof]
        self.dda = acc
        self.F = self.forces()
        
        self._updated = True
        
    def forces(self):
        """This function returns the forces for Newmark
        """
        C = self.control
        F = np.zeros(self.dof) # modal force
        ones = np.ones(self._mesh.size)
        dforce = np.zeros(self._mesh.size) # Distributed forces
        
        # Sum all the distributed forces
        for f in C['forces']:
            if f['type'] == 'gravity':
                dforce += C['m'] * f['value'] * ones
            elif f['type'] == 'distributed':
                obj = self._mesh.OREG.get(f['object'])  # Get the object from the registry
                try:
                    dforce += obj.beamDistributedForce(self)  # The object must have a beamForce method
                except:
                    print("---> ERROR: Object " + obj.name + " does not have a beamDistributedForce method.")
            if f['type'] == "boundary":
                for boundary in self.boundary:
                    dforce += boundary.f
                    
        # Assemble the distributed forces
        for i in range(self.dof):
            F[i] += si.simps(dforce * self.eigen.vectors[i], self._mesh.x)
            
                #  Add the concentrated forces
        for f in C['forces']:
            if f['type'] == 'concentrated':
                obj = self._mesh.OREG.get(f['object'])
                loc = f['location']
                pforce = obj.beamConcentratedForce(self)
                for i in range(self.dof):
                    evector = self.eigen.vectors[i]
                    F[i] += pforce * evector[loc] # Dirac Delta Integral

        return F


    def RHS(self, time, state):
        Q = np.zeros(self.sof)
        # State Space Force
        Minv = self.M  # Inverse of the mass matrix
        Q[self.dof:self.sof] = np.dot(Minv, self.forces())
        # Calculate the RHS
        self.rhs = np.dot(self.S, state) + Q # Store the acceleration

        return self.rhs

    # ----- Private methods ----- #
    def _createModalFunctions(self):
        # Mass function
        self.m = self.control['m'] * self.eigen.vectors
        # Damping function using modal damping
        if self.control['solution']['damping']['type'] == 'modal':
            for i in range(self.dof):
                self.c[i] = self.control['m'] * (
                            self.control['solution']['damping']['ratios'][i] *
                            self.eigen.values[i] * self.eigen.vectors[i])
        else:
            print("---> WARNING! Damping type not found. No damping applied")

        # Stiffness function
        self.k = self.control['k']

    def _createStateSpaceMatrix(self):
        # State matrix assumming mass normalized eigenvectors
        dof, sof = self.dof, self.sof
        I = np.identity(dof)
        self.M = I  # Linear mass
        self.Minv = self.M
        
        # Assemble the linear stiffness for self-adjoint operator
        for i in range(dof):
            self.K[i, i] = self.eigen.values[i] ** 2

        # Assemble the damping matrix
        if self.control['solution']['damping']['type'] == 'Rayleigh':
            print("---> Rayleligh damping cannot be used with only one mode. No damping applied")
        elif self.control['solution']['damping']['type'] == 'modal':
            for i in range(dof):
                self.C[i, i] = 2 * self.control['solution']['damping']['ratios'][i] * self.eigen.values[i]
        else:
                print("---> WARNING!. Damping type not found. No dampign applied")

        # State matrix
        self.S[0:dof, dof:sof] = I
        self.S[dof:sof, 0:dof] = -np.dot(self.Minv, self.K)
        self.S[dof:sof, dof:sof] = -np.dot(self.Minv, self.C)

    def _createEigenSystem(self, control):
        # Create the eigen system object
        eigenValues = self._eigenValues(control)
        eigenVectors = self._eigenVectors(eigenValues)
        self.eigen = esys.eigenSystemVector(eigenValues, 
                                            eigenVectors,
                                            calculate=True)
        freqs = [f/6.2832 for f in self.freqs()]
        printArray(freqs, name="System frequencies (Hz)", file="./solid/beamFrequencies.out")
        # [plt.plot(self._mesh.x, ev) for ev in eigenVectors] # Plot Modes
        plt.show()



    # Getters
    def freqs(self):
        return self.eigen.values

    def mesh(self):
        return self._mesh

    def control(self):
        return self.control

    def _messages(self, n):
        pass

    def _eigenValues(self, control):
        # Create the eigenvalue list (based on the frequency)
        eigenValues = []
        k, m = self.control['k'], self.control['m']
        eigenValues.append(eval.eigenValue(k / m))

        return eigenValues

    def _eigenVectors(self, eigenValues):

        eigenVectors = []
        eigenVectors.append(evec.eigenVector(np.array(1), 
                                             mesh=self._mesh, 
                                             calculate=True, 
                                             normalize='mass')
                            )
        
        return eigenVectors
