# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#    p    #     version: 0.2
#    y    #     date: 03/11/2023
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   The class  works either as an analytical model or with a Calculix model
#   The model type is selected from the key "method" in the solid-solution dict
#   Supports reconstruction of displacements
#   Modal and Rayleigh damping models are available
#   Both 1st and 2nd order time integrator can be used with this beam
# Warnings:
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

class flextensionalBeam(solidModel):
    def __repr__(self):
        return 'beamModel'
    
    def __init__(self, execution, control, mesh, boundary, time):
        super().__init__(execution, control, mesh, boundary, time)  # Call the base class
        # ----- Public attributes ----- #
        # Degrees of freedom
        self.dof = control['solution']['modes']     # Modal DOFs (one for the base motion)
        self.sof = 2 * self.dof                         # State Space DOFS
        self.eigen = None                               # Eigensystem
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
            print("     WARNING: Beam initial conditions set to zero. ")

        # ----- Procedures ----- #
        # Add calculated parameters to the control dictionary
        control['L'] = self._mesh.L  # Add the length to the parameters dictionary
        control['I'] = (control['section']['b'] * (control['section']['h'] ** 3) / 12)
        control['A'] = control['section']['b'] * control['section']['h']
        control['rho'] = control['material']['rho'] 
        control['m'] = control['A'] * control['rho']  # mass per unit
        control['EI'] = control['material']['E'] * control['I']
        control['hm'] = 0.5 * self.control['section']['h']  # Mid-distance

        # Create the eigensystem (it doet not use boundary info, just the mesh)
        self._createEigenSystem(control)

        # Call forces only if magnetic stiffness is set up
        for f in control['forces']:
            if f['type'] == 'distributedLineal':
                self.forces(True)

                freqs = [f/6.2832 for f in self.freqs()]
                printArray(freqs, name="Beam frequencies modified by magnetic stiffness(Hz)",
                           file="./solid/beamMagneticFrequencies.out")
                plt.show()
        
        # Calculate the dimensionless numbers
        self.calcNumbers()
        
        # Calculate the state matrix S
        self._createMatrices()  
        
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
        self.update(None, initialState, self.dda, boundary=False)
        

    # Update the kinematics without updating the state (no convergence)
    def update(self, time, state, acc, boundary=True):
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
        if boundary:
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
             
        
        # Elastic, kinetic and potential energies (respectively)
        # self.energy = np.zeros(4)
        # grav = 0.0
        # for i in C['forces']:
        #     if i['type'] == 'gravity':
        #         grav = i['value']
        #         break
        # self.energy[0] = 0.5 * si.simps(C['EI'] * ddx**2, x)
        # self.energy[1] = 0.5 * si.simps(C['rho'] * C['A'] * self.dy['mid']**2, x)
        # self.energy[2] = -si.simps(C["rho"] * grav * C['A'] * self.v['mid'], x)
        # # Add th boundary condition energies
        # if ((C['bc']['type'] == 'elastic-free') or 
        #     (C['bc']['type'] == 'spring-free')):
        #     self.energy[0] += 0.5 * C['bc']['Kt'] * self.v['mid'][0]**2
        #     self.energy[1] += 0.5 * C['bc']['Mt'] * self.dy['mid'][0]**2
        #     self.energy[2] -= C['bc']['Mt'] * grav * self.v['mid'][0]
        # # Total energy and power
        # self.energy[3] = self.energy[0] + self.energy[1] + self.energy[2]
        # self.power = (self.energy - self.energy0) / self._time.delta
        # # Store the energy
        # self.energy0 = self.energy.copy()
        
        self._updated = True
        
    def forces(self, magnetStiffness=False):
        """This function returns the forces for Newmark
        """
        C = self.control
        F = np.zeros(self.dof) # modal force
        ones = np.ones(self._mesh.size)
        dforce = np.zeros(self._mesh.size) # Distributed forces
        
        # Sum all the distributed forces
        if magnetStiffness == True:
            for f in C['forces']:
                if f['type'] == 'distributedLineal':
                    # Get the object from the registry
                    obj = self._mesh.OREG.get(f['object'])
                    try:
                        # The object must have a beamForce method
                        dforce += obj.beamDistributedForce(self)
                    except:
                        print("---> ERROR: Object " + obj.name +
                              " does not have a beamDistributedForce method.")

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
                loc = f['location'] # Degree of freedom where the force is located (node index)
                if obj.type == 'spring':
                    pforce = obj.force(self.v['mid'][f['location']])
                if obj.type == 'damper':
                    pforce = obj.force(self.dy['mid'][f['location']])
                if obj.type == 'beamConcentratedForce':
                    pforce = obj.beamConcentratedForce(self)
                for i in range(self.dof):
                    evector = self.eigen.vectors[i]
                    F[i] += pforce * evector[loc] # Dirac Delta Integral
        # Add the force for the point mass in the BC (only if gravity is present)
        if ((C['bc']['type'] == 'elastic-free') or 
            (C['bc']['type'] == 'spring-free')):
            for f in C['forces']:
                if f['type'] == 'gravity':
                    bcforce = C['bc']['Mt'] * f['value']     
                    for i in range(self.dof):
                        # Dirac Delta Integral
                        F[i] += bcforce * self.eigen.vectors[i][0] 

        return F


    def RHS(self, time, state):
        Q = np.zeros(self.sof)
        # State Space Force
        Minv = self.M  # Inverse of the mass matrix
        Q[self.dof:self.sof] = np.dot(Minv, self.forces())
        # Calculate the RHS
        self.rhs = np.dot(self.S, state) + Q # Store the acceleration

        return self.rhs
    

    def calcNumbers(self):
        # Fill the reference parameters
        self.lRef = self.control['L']
        self.tRef = self.control['section']['h']
        self.uRef = self.control['L'] / 10.0  # Approximate reference displacement as L/10, no justification for this!
        # We take the largest natural frequency to calculate the reference velocity
        # Scaling with the sound speed in the solid would be not very realistic
        self.vRef = self.freqs()[-1] * self.uRef

        # Calculate the dimensionless numbers
        super().calcNumbers()

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
        self.k = self.control['material']['E'] * self.control['I'] * self.eigen.d4

    def _createMatrices(self):
        # State matrix assumming mass normalized eigenvectors
        dof, sof = self.dof, self.sof
        I = np.identity(dof)
        self.M = I  # Linear mass
        self.Minv = self.M

        # Assemble the linear stiffness for the self-adjoint operator
        for i in range(dof):
            self.K[i, i] = self.eigen.values[i] ** 2

        # Assemble the damping matrix
        if self.control['solution']['damping']['type'] == 'Rayleigh':
            print("---> Rayleigh damping was chosen...")
            if dof > 1:
                s1 = self.control['solution']['damping']['s1']
                s2 = self.control['solution']['damping']['s2']
                w1 = self.eigen.values[0]
                w2 = self.eigen.values[1]
                alpha = 2 * ( s2 * w1**2 * w2 - s1 * w1 * w2**2) / (w1**2 - w2**2)
                beta = 2 * ( s1 * w1 - s2 * w2) / (w1**2 - w2**2)
                for i in range(dof):
                    self.C[i, i] = alpha * self.M[i, i] + beta * self.K[i, i]
            else:
                print("---> Rayleligh damping cannot be used with only one mode. No damping applied")
        elif self.control['solution']['damping']['type'] == 'modal':
            print("---> Modal damping was chosen...")
            for i in range(dof):
                self.C[i, i] = 2 * self.control['solution']['damping']['ratios'][i] * self.eigen.values[i]
            printArray(self.C, name="Damping Matrix", file="./solid/"+self.name+"_dampingMatrix.out")
        else:
                print("---> WARNING!. Damping type not found. No dampign applied")


        # State Matrix
        self.S[0:dof, dof:sof] = I
        self.S[dof:sof, 0:dof] = -np.dot(self.Minv, self.K)
        self.S[dof:sof, dof:sof] = -np.dot(self.Minv, self.C)

    def _createEigenSystem(self, control):
        # Create the eigen system object
        if control['solution']['type'] == 'modal':
            if control['solution']['method'] == 'analytic':
                eigenValues = self._eigenValues(control)
                eigenVectors = self._eigenVectors(eigenValues)
                self.eigen = esys.eigenSystemVector(eigenValues, 
                                                    eigenVectors,
                                                    calculate=True)
                freqs = [f/6.2832 for f in self.freqs()]
                printArray(freqs, name="Beam frequencies (Hz)", file="./solid/"+self.name+"_beamFrequencies.out")
                # [plt.plot(self._mesh.x, ev) for ev in eigenVectors] # Plot Modes
                plt.show()
            elif control['solution']['method'] == 'calculix':
                import pyFSI.models.solidModels.calculixBeam as cx
                self.eigen = self._calculixEigenSystem(control['solution']['modes'])
            else:
                print("--> ERROR: Beam solution method " + control['solution']['method'] + " not kwnown !")


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
        nBetas, bc = control['solution']['modes'], control['bc']
        L, EI, m = control['L'], control['EI'], control['m']

        if bc['type'] == "clamped-free":
            betaL = np.array([1.875104068711961, 4.694091132974175, 7.854757438237612,
                              10.995540734875467, 14.13716839104647, 17.278759532088234])
            betas = betaL[0:nBetas] / L

        elif bc['type'] == "supported-supported":
            betaL = np.array([3.141592653589793, 6.283185307179586, 9.42477796076938,
                              12.566370614359172, 15.707963267948966, 18.84955592153876])
            betas = betaL[0:nBetas] / L

        elif bc['type'] == "clamped-supported":
            betaL = np.array([3.9265975952148438, 7.068580627441406, 10.210182189941406,
                              13.351768493652344, 14.137166976928711, 16.49335479736328])
            betas = betaL[0:nBetas] / L

        elif bc['type'] == "spring-free": 
            with warnings.catch_warnings():  # Ignore warings for this part of the code (overflow is normal)
                warnings.simplefilter("default")
                Kt, Ct, Mt = bc['Kt'], bc['Ct'], bc['Mt']
                # beta function RHS
                # Method using my routine (insert one bc into the other)
                fvalue0 = lambda b: ((-2*(-(Kt*m) + b**4*EI*Mt + (-(Kt*m) + b**4*EI*Mt)*(1/np.cos(b*L)) * (1/np.cosh(b*L)) 
                                    +  b**3*EI*m*(np.tan(b*L) + np.tanh(b*L))))/(m*(np.tan(b*L) + np.tanh(b*L))))
                # Method using the nullspace (Santi's method). Currently using it.
                fvalue1 = lambda b: ((2*b**6*(-(Kt*m) + EI*Mt*b**4 + np.cosh(L*b)*((-(Kt*m) + EI*Mt*b**4)*np.cos(L*b) 
                                     + EI*m*b**3*np.sin(L*b)) + EI*m*b**3*np.cos(L*b)*np.sinh(L*b)))/m)
                                    
                fprime0 =  lambda b: ( (-2*(4*b**3*EI*Mt + 4*b**3*EI*Mt*(1/(np.cos(b*L)))*(1/(np.cosh(b*L))) + 
                                    b**3*EI*m*(L*(1/(np.cos(b*L)))**2 + L*(1/(np.cosh(b*L)))**2) + 
                                    L*(-(Kt*m) + b**4*EI*Mt)*(1/(np.cos(b*L)))*(1/(np.cosh(b*L)))*np.tan(b*L) - 
                                    L*(-(Kt*m) + b**4*EI*Mt)*(1/(np.cos(b*L)))*(1/(np.cosh(b*L)))*np.tanh(b*L) + 
                                    3*b**2*EI*m*(np.tan(b*L) + np.tanh(b*L))))/(m*(np.tan(b*L) + np.tanh(b*L))) + 
                                    (2*(L*(1/(np.cos(b*L)))**2 + L*(1/(np.cosh(b*L)))**2)*
                                    (-(Kt*m) + b**4*EI*Mt + (-(Kt*m) + b**4*EI*Mt)*(1/(np.cos(b*L)))*(1/(np.cosh(b*L))) + 
                                    b**3*EI*m*(np.tan(b*L) + np.tanh(b*L))))/(m*(np.tan(b*L) + np.tanh(b*L))**2))
                
                # Method using a function to make it work with Newton's method (requires a list as output)
                def fvalue(b):
                    # Beta function
                    f = ((-2*(-(Kt*m) + b**4*EI*Mt + (-(Kt*m) + b**4*EI*Mt)*(1/np.cos(b*L)) * (1/np.cosh(b*L)) 
                                        +  b**3*EI*m*(np.tan(b*L) + np.tanh(b*L))))/(m*(np.tan(b*L) + np.tanh(b*L))))
                    # First derivative
                    d = ( (-2*(4*b**3*EI*Mt + 4*b**3*EI*Mt*(1/(np.cos(b*L)))*(1/(np.cosh(b*L))) + 
                                        b**3*EI*m*(L*(1/(np.cos(b*L)))**2 + L*(1/(np.cosh(b*L)))**2) + 
                                        L*(-(Kt*m) + b**4*EI*Mt)*(1/(np.cos(b*L)))*(1/(np.cosh(b*L)))*np.tan(b*L) - 
                                        L*(-(Kt*m) + b**4*EI*Mt)*(1/(np.cos(b*L)))*(1/(np.cosh(b*L)))*np.tanh(b*L) + 
                                        3*b**2*EI*m*(np.tan(b*L) + np.tanh(b*L))))/(m*(np.tan(b*L) + np.tanh(b*L))) + 
                                        (2*(L*(1/(np.cos(b*L)))**2 + L*(1/(np.cosh(b*L)))**2)*
                                        (-(Kt*m) + b**4*EI*Mt + (-(Kt*m) + b**4*EI*Mt)*(1/(np.cos(b*L)))*(1/(np.cosh(b*L))) + 
                                        b**3*EI*m*(np.tan(b*L) + np.tanh(b*L))))/(m*(np.tan(b*L) + np.tanh(b*L))**2))
                    return [f, d]
                betas = []
                nroot, step, betai, betaTol, betaMax = 0, 2, 0.1, 1E-3, 5E3  # Root finding controls
                while nroot < nBetas:
                    betaf = betai + step
                    if np.sign(fvalue1(betai)) != np.sign(fvalue1(betaf)):
                        sol = optimize.root_scalar(fvalue1, bracket=[betai, betaf], rtol=1E-15, method='brentq', fprime=False)
                        if sol.root > (betai + betaTol * betai): # Test for repeated roots (also avoids the root beta=0)
                            #root = optimize.root_scalar(fvalue, x0=sol.root, rtol=1E-15, fprime=True,  method='newton')
                            betas.append(sol.root)
                           # betas.append(sol.root) # Store the root found
                            nroot += 1
                    betai += step
                    if betai > betaMax:
                        sys.exit("--> ERROR: Maximum beta reached. Current roots are: " + str(betas))
            warnings.simplefilter("default")

        elif bc['type'] == "elastic-free":
            with warnings.catch_warnings():  # Ignore warings for this part of the code (overflow is normal)
                warnings.simplefilter("ignore")
                Kt, Ct, Mt = bc['Kt'], bc['Ct'], bc['Mt']
                Kr, Cr, Mr = bc['Kr'], bc['Cr'], bc['Mr']
                EIm = EI * m
                # beta function RHS using the nullspace method
                fvalue = lambda b: ((2*b**6/m**2)*(-(Kr*Kt*m**2) + EIm*(EIm + Kt*Mr + Kr*Mt)*b**4 - EI**2*Mr*Mt*b**8 +
                                    np.cosh(L*b)*(-((Kr*Kt*m**2 + EIm*(EIm - Kt*Mr - Kr*Mt)*b**4 +
                                    EI**2*Mr*Mt*b**8)*np.cos(L*b)) + EIm*b*(-(Kt*m) + Kr*m*b**2 + EI*Mt*b**4 -
                                    EI*Mr*b**6)*np.sin(L*b)) + EIm*b*(Kt*m + Kr*m*b**2 -
                                    EI*b**4*(Mt + Mr*b**2))*np.cos(L*b)*np.sinh(L*b)))

                # b derivative. It does not work with root finding, I do not know why (error in the copy from Mathematica?)
                fprime = lambda b:  ((2*b**5*(-6*Kr*Kt*m**2 + 10*EI*m*(Mr*Kt + EI*m + Kr*Mt)*b**4 - 14*EI**2*Mr*Mt*b**8 +
                                    np.cosh(L*b)*(-2*(3*Kr*Kt*m**2 + EI*m*(-5*Mr*Kt + 5*EI*m - Kr*L*m - 5*Kr*Mt)*b**4 +
                                    EI**2*Mr*(L*m + 7*Mt)*b**8)*np.cos(L*b) + b*(Kr*Kt*L*m**2 + EI**2*b**4*(L*m**2 + 11*m*Mt -
                                    13*Mr*m*b**2 + Mr*L*Mt*b**4) - EI*m*(7*Kt*m - 9*Kr*m*b**2 + L*(Mr*Kt + Kr*Mt)*b**4)) *
                                    np.sin(L*b)) + b*(-((Kr*Kt*L*m**2 + EI**2*b**4*(L*m**2 + 11*m*Mt + 13*Mr*m*b**2 +
                                    Mr*L*Mt*b**4) - EI*m*(7*Kt*m + 9*Kr*m*b**2 + L*(Mr*Kt + Kr*Mt)*b**4))*np.cos(L*b)) +
                                    2*EI*L*m*b*(-(Kt*m) + EI*Mt*b**4)*np.sin(L*b))*np.sinh(L*b)))/m**2)

                betas = []
                nroot, step, betai, betaTol, betaMax = 0, 2, 0.1, 1E-3, 5E3  # Root finding controls
                while nroot < nBetas:
                    betaf = betai + step
                    if np.sign(fvalue(betai)) != np.sign(fvalue(betaf)):
                        sol = optimize.root_scalar(fvalue, bracket=[betai, betaf], method='brentq', fprime=False)
                        if sol.root > (betai + betaTol * betai): # Test for repeated roots (also avoids the root beta=0)
                            betas.append(sol.root) # Store the root found
                            nroot += 1
                        # Newton Raphson loop (not needed)
                        # root = optimize.root_scalar(fvalue, x0=sol.root, fprime=True, rtol=1E-15,  method='newton')
                        # roots.append(root.root)
                    betai += step
                    if betai > betaMax:
                        sys.exit("--> ERROR: Maximum beta reached. Current roots are: " + str(betas))
                warnings.simplefilter("default")

        else:
            bcTypes = ["clamped-free", "supported-supported", "clamped-supported", "elastic-free", "spring-free"]
            print("--> ERROR: Non-recognized beam BC: " + str(control['bc']['type']))
            print(    "Available types are: ", bcTypes)
            sys.exit()
        print("     INFO: The beta equation has converged. Betas found are: ")
        #printArray(betas, name="Betas")
        
        self.betas = betas

        # Create the eigenvalue list (based on the frequency)
        eigenValues = []
        for j in self.betas:
            eigenValues.append(eval.eigenValue((EI/m)**0.5 * j**2)) # Create the eigenvalue (omega) from beta

        return eigenValues

    def _eigenVectors(self, eigenValues):
        control, bc = self.control, self.control['bc']
        x, L, EI, m = self._mesh.x, control['L'], control['EI'], control['m']
        eigenVectors = []

        for b in self.betas:
            bx, bL, bEIm = b * x,  b * L, b * EI * m
            # Select the type of bc
            if self.control['bc']['type'] == "clamped-free":
                Mt = 0
                vector = (np.cosh(bx) - np.cos(bx) - ((np.sinh(bL) - np.sin(bL)) / (np.cosh(bL) + np.cos(bL)))
                          * (np.sinh(bx) - np.sin(bx)))
                eigenVectors.append(evec.eigenVector(vector, 
                                                     mesh=self._mesh,
                                                     calculate=True,
                                                     normalize='mass', 
                                                     mass=m))

            elif self.control['bc']['type'] == "supported-supported":
                normalization = 'mass'
                Mt = 0
                vector = np.sin(bx)
                eigenVectors.append(evec.eigenVector(vector, 
                                                     mesh=self._mesh,
                                                     calculate=True,
                                                     normalize='mass', 
                                                     mass=m))


            elif self.control['bc']['type'] == "clamped-supported":
                normalization = 'mass'
                Mt = 0
                vector = ((np.cos(bx) - np.cosh(bx)) - ((np.cos(bL) - np.cosh(bL)) / (np.sin(bL) - np.sinh(bL)))
                          * (np.sin(bL) - np.sinh(bL)))
                eigenVectors.append(evec.eigenVector(vector, 
                                                     mesh=self._mesh,
                                                     calculate=True,
                                                     normalize='mass', 
                                                     mass=m))


            elif self.control['bc']['type'] == "spring-free":
                Kt, Ct, Mt = bc['Kt'], bc['Ct'], bc['Mt']
                sinbx, sinhbx, sinbL, sinhbL = np.sin(bx), np.sinh(bx), np.sin(bL), np.sinh(bL)
                cosbx, coshbx, cosbL, coshbL = np.cos(bx), np.cosh(bx), np.cos(bL), np.cosh(bL)
                tanbx, tanhbx, tanbL, tanhbL = np.tan(bx), np.tanh(bx), np.tan(bL), np.tanh(bL)
                secbL, sechbL = 1/(cosbL),  1/coshbL 
                
                # Mode (non-normalized)
                vector = (np.sin(bx) - np.sinh(bx) + (np.cos(bx)*(1 + 
                secbL*sechbL - np.tan(bL)*np.tanh(bL)))/(np.tan(bL) + 
                np.tanh(bL)) + (np.cosh(bx)*(1 + secbL*sechbL + np.tan(bL) * 
                np.tanh(bL)))/(np.tan(bL) + np.tanh(bL)))
                
                # First derivative of the mode (non-normalized)
                first = (b*(-coshbx + (-sinbx + sinhbx + secbL*sechbL*(-sinbx + 
                sinhbx) + cosbx*tanbL + (cosbx + (sinbx + sinhbx)*tanbL)*tanhbL)/
                (tanbL + tanhbL)))

                second = ((b**2*(-((cosbx - coshbx)*(1 + secbL*sechbL)) - 
                (sinbx + sinhbx)*tanbL - (sinbx + sinhbx - (cosbx + coshbx) * 
                tanbL)*tanhbL))/(tanbL + tanhbL))

                third =  ((b**3*((1 + secbL*sechbL)*(sinbx + sinhbx) - 
                (cosbx + coshbx)*tanbL - (cosbx + coshbx + (sinbx - sinhbx) * 
                tanbL)*tanhbL))/(tanbL + tanhbL))

                fourth = ((b**4*((cosbx + coshbx)*(1 + secbL*sechbL) + 
                (sinbx - sinhbx)*tanbL + (sinbx - sinhbx + (-cosbx + 
                coshbx)*tanbL)*tanhbL))/(tanbL + tanhbL))

                vectorIx = (-(cosbx/b) - coshbx/b + (sinbx*(1 + secbL*sechbL 
                - tanbL*tanhbL))/(b*(tanbL + tanhbL)) + (sinhbx*(1 + secbL * 
                sechbL + tanbL*tanhbL))/(b*(tanbL + tanhbL)))

                vectorIL = 2 / b
                
                # Calculation of the normalization factor
                deno = si.simps(m * vector**2, x) + Mt * vector[0]**2
                if deno < 0:  # We need to invert the eigenvector to get a positive eigenvalue
                    vector *= -1
                    deno *= -1
                    print("    WARNING: Mirroring vector with beta ", b)
                normFactor = (1 / deno)**0.5

                # Create the eigenvector object (normalization is done inside)
                eigenVectors.append(evec.eigenVector(vector, 
                                    mesh=self._mesh,
                                    calculate=True,
                                    normalize='general', 
                                    factor=normFactor,
                                    first=first,
                                    second=second,
                                    third=third,
                                    fourth=fourth,
                                    integralx=vectorIx,
                                    integralL=vectorIL))

            elif self.control['bc']['type'] == "elastic-free":
                Kt, Ct, Mt = bc['Kt'], bc['Ct'], bc['Mt']
                Kr, Cr, Mr = bc['Kr'], bc['Cr'], bc['Mr']
                sinbx, sinhbx, sinbL, sinhbL = np.sin(bx), np.sinh(bx), np.sin(bL), np.sinh(bL)
                cosbx, coshbx, cosbL, coshbL = np.cos(bx), np.cosh(bx), np.cos(bL), np.cosh(bL)
                tanbx, tanhbx, tanbL, tanhbL = np.tan(bx), np.tanh(bx), np.tan(bL), np.tanh(bL)

                vector = (sinbx + (cosbx*(b**4*EI*Mr - Kr*m + coshbL*((b**4*EI*Mr - Kr*m)*cosbL - bEIm*sinbL) +
                        (bEIm*cosbL + (-(b**4*EI*Mr) + Kr*m)*sinbL)*sinhbL))/
                        (-bEIm + coshbL * (bEIm * cosbL + (b ** 4 * EI * Mr - Kr * m) * sinbL) +
                         ((b**4*EI*Mr - Kr*m)*cosbL + bEIm*sinbL) * sinhbL) -
                        ((-bEIm + coshbL*(bEIm*cosbL + (b**4*EI*Mr - Kr*m)*sinbL) +
                        ((b**4*EI*Mr - Kr*m)*cosbL - bEIm*sinbL)*sinhbL)*sinhbx)/
                        (-bEIm + coshbL * (bEIm * cosbL + (b ** 4 * EI * Mr - Kr * m) * sinbL) +
                         ((b**4*EI*Mr - Kr*m)*cosbL + bEIm*sinbL) * sinhbL) +
                        (np.cosh(bx)*(b**4*EI*Mr*cosbL - Kr*m*cosbL + (b**4*EI*Mr - Kr*m)*(1/coshbL) - bEIm*sinbL +
                        (bEIm*cosbL + (b**4*EI*Mr - Kr*m)*sinbL)*tanhbL))/
                        (bEIm*cosbL - bEIm*(1/coshbL) + b**4*EI*Mr*sinbL - Kr*m*sinbL +
                        ((b**4*EI*Mr - Kr*m)*cosbL + bEIm*sinbL)*tanhbL))

        #         SechbL = (1/np.cosh(bL))
        #         vector0 =  (np.sin(b*x) + (np.cos(b*x)*(-(Kr*m) + b*4*EI*Mr - 
        #       np.cosh(b*L)*((Kr*m - b*4*EI*Mr)*np.cos(b*L) + b*EI*m*np.sin(b*L)) + 
        #       (b*EI*m*np.cos(b*L) + (Kr*m - b**4*EI*Mr)*np.sin(b*L))*np.sinh(b*L)))/
        #   (-(b*EI*m) + np.cosh(b*L)*
        #      (b*EI*m*np.cos(b*L) + (-(Kr*m) + b**4*EI*Mr)*np.sin(b*L)) + 
        #     ((-(Kr*m) + b**4*EI*Mr)*np.cos(b*L) + b*EI*m*np.sin(b*L))*np.sinh(b*L)) - 
        #  ((-(b*EI*m) + np.cosh(b*L)*
        #        (b*EI*m*np.cos(b*L) + (-(Kr*m) + b**4*EI*Mr)*np.sin(b*L)) - 
        #       ((Kr*m - b*4*EI*Mr)*np.cos(b*L) + b*EI*m*np.sin(b*L))*np.sinh(b*L)) *
        #     np.sinh(b*x))/
        #   (-(b*EI*m) + np.cosh(b*L)*
        #      (b*EI*m*np.cos(b*L) + (-(Kr*m) + b**4*EI*Mr)*np.sin(b*L)) + 
        #     ((-(Kr*m) + b**4*EI*Mr)*np.cos(b*L) + b*EI*m*np.sin(b*L))*np.sinh(b*L)) + 
        #  (np.cosh(b*x)*(-(Kr*m*np.cos(b*L)) + b*4*EI*Mr*np.cos(b*L) + 
        #       (-(Kr*m) + b**4*EI*Mr)*SechbL - b*EI*m*np.sin(b*L) + 
        #       (b*EI*m*np.cos(b*L) + (-(Kr*m) + b**4*EI*Mr)*np.sin(b*L))*np.tanh(b*L)))/
        #   (b*EI*m*np.cos(b*L) - b*EI*m*SechbL - Kr*m*np.sin(b*L) + 
        #     b**4*EI*Mr*np.sin(b*L) + 
        #     ((-(Kr*m) + b**4*EI*Mr)*np.cos(b*L) + b*EI*m*np.sin(b*L))*np.tanh(b*L)))

                dvector = (b*(cosbx - (np.cosh(bx)*(-bEIm + coshbL*(bEIm*cosbL +
                       (-(Kr*m) + b**4*EI*Mr)*sinbL) - ((Kr*m- b**4*EI*Mr)*cosbL + bEIm*sinbL)*sinhbL)) /
                       (-bEIm + coshbL*(bEIm*cosbL + (-(Kr*m) + b**4*EI*Mr)*sinbL) +
                       ((-(Kr*m) + b**4*EI*Mr)*cosbL + bEIm*sinbL)*sinhbL) +
                       (sinbx*(Kr*m - b**4*EI*Mr + coshbL*((Kr*m - b**4*EI*Mr)*cosbL + bEIm*sinbL) +
                       (-(bEIm*cosbL) + (-(Kr*m) + b**4*EI*Mr)*sinbL)*sinhbL))/(-bEIm + coshbL *
                       (bEIm*cosbL + (-(Kr*m) + b**4*EI*Mr)*sinbL) +
                       ((-(Kr*m) + b**4*EI*Mr)*cosbL + bEIm*sinbL)*sinhbL) +
                       (sinhbx*(-(Kr*m*cosbL) + b**4*EI*Mr*cosbL +
                       (-(Kr*m) + b**4*EI*Mr)*(1/coshbL) - bEIm*sinbL +
                       (bEIm*cosbL + (-(Kr*m) + b**4*EI*Mr)*sinbL)*tanhbL)) /
                       (bEIm*cosbL - bEIm*(1/coshbL) - Kr*m*sinbL +
                       b**4*EI*Mr*sinbL + ((-(Kr*m) + b**4*EI*Mr)*cosbL + bEIm*sinbL)*tanhbL)))
                # append the normalized eigenvector to the list
                deno = si.simps(m * vector**2, x) + Mt * vector[0]**2 - Mr * dvector[0]**2
                # if deno < 0:  # We need to invert the eigenvector to get a positive eigenvalue
                #     vector *= -1
                #     deno *= -1
                #     print("    INFO: Mirroring vector with eigenvalue ", ev)
                normFactor = (1 / deno)**0.5
                eigenVectors.append(evec.eigenVector(vector, 
                                                     mesh=self._mesh, 
                                                     calculate=True, 
                                                     normalize='general', 
                                                     factor=normFactor))

        J = np.zeros((self.dof, self.dof))
        for i in range(self.dof):
            for j in range(self.dof):
                J[i, j] = (si.simps(m*eigenVectors[i]*eigenVectors[j], self._mesh.x) 
                          + Mt * eigenVectors[i][0]*eigenVectors[j][0])
        
        printArray(J, name="Orthogonality", file="./solid/"+self.name+"_orthogonality.out")
        
        return eigenVectors
