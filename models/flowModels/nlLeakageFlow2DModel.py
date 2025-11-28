# --------------------------------------------------------------------------- #
#    p    #     version: 0.2
#    y    #     date: 11/03/2024
#    F    #     author: Martin Saravia
#    S    #     description: NL Leakage Flow Class
#    I    #     return: flow object
# --------------------------------------------------------------------------- #
# Notes:
#   Nonlinear leakage flow model.
#
# New version:
#   Adds the initialize method to avoid recalculation of the W operators that
#   are a function of the solid. This implies a new RHS method. The method is
#   called by the solver at the beginning of each time step.
#
# Warnings:
#   Local Re is the same for all regions. Not as above Eq. 2.30
#   eta calculation must be checked
#   Check the consideration of the beam thickness (hardcoded to 1)
#   Only two regions allowed
#   Regions has just one flexible boundary
#   Inconsistency: self.control['bc']['inlet']['zeta'] is the same for regions
# --------------------------------------------------------------------------- #

from abc import ABC
import numpy as np
import scipy.integrate as si

from pyFSI.mesh.region.fsiRegion1D import fsiRegion1D
from pyFSI.models.flowModels.flowBase import flowModel
from pyFSI.models.properties.boundaryLayer import boundaryLayer
from pyFSI.vectors.eigen.eigenVector import eigenVector
from pyFSI.models.properties.dimensionlessNumbers import ReynoldsNumber
from pyFSI.fields.boundary import boundaryConditions
from pyFSI.post.prints import *
from threading import Thread

class nlLeakageFlow2D(flowModel, ABC):
    def __repr__(self):
        return 'nonlinear leakageFlow2DModel '

    def __init__(self, execution, control, mesh, boundary, time, name='NN'):
        super().__init__(execution, control, mesh, boundary, time)

        # ----- Public attributes ----- #
        self.dof = len(self.regions)                # Number of Q equations (1 per region)
        self.Q0 = np.zeros(self.dof)                # Flow rate at x = 0
        self.vx = np.zeros(self.dof, dtype=object)
        self.dQ0 = np.zeros(self.dof) # 
        self.Q = np.zeros(self.dof, dtype=object)
        self.dQ = np.zeros(self.dof, dtype=object)  # Flow rate time derivative
        self.px = np.zeros(self.dof, dtype=object)  # Pressure distribution
        self.deltaPx = np.zeros(1, dtype=object)
        self.totalForce = 0
        self.Forces = np.zeros((mesh.size*2, 3))    # The 3D Forces
        self.Dp = None                              # Pressure difference: pOut - pIn
        self.power = {
            'inertial': 0.0,
            'convective': 0.0,
            'convective2': 0.0,
            'pressure': 0.0,
            'viscous': 0.0, 
            'heat': 0.0
        }
        self.converged = [False] * self.dof
        # State variables
        self.state = self.Q0

        # Output variables mapping
        self.varMap['flowRates'] = 'Q0'
        self.varMap['pressures'] = 'px'
        self.varMap['flowSpeeds'] = 'vRef'
        self.varMap['inletVelocity'] = 'vi'
        self.varMap['power'] = 'power'

        # ----- Private attributes ----- #
        self._ti = None  # Current time
        self._th = 1  # Thickness of the model
        self._pIn = None  # Inlet pressure Object
        self._pOut = None  # Outlet pressure
        self._zetaIn = None  # Inlet loss factor
        self._zetaOut = None  # Outlet loss factor
        self._pTol = 1E-4  # Pressure convergence tolerance
        # self._f0 = np.zeros(self.dof)  # Viscous friction factor
        # self._xix = np.zeros(self.dof)  # Nonlinear profile factor
        # self._eta = np.zeros(self.dof)  # Derivative of f(Qx) at qx0
        self.dRef = np.zeros(self.dof)  # Reference width
        self.lRef = mesh.L  # Fixed Reference length
        self.vRef = np.zeros(self.dof)  # Reference velocity
        self.eRef = np.zeros(self.dof)  # Reference something
        self.f0 = np.zeros(self.dof)
        self.xix = np.zeros(self.dof)
        self.eta = np.zeros(self.dof)
        self.Rd = [None] * self.dof  # Emtpy list of Reynolds objects
        # Size of the associated eigensystem
        # self._gDof = self.regions[0].eigen().size

        # Initialize the boundary conditions
        # Pressure BC
        pBCDict = self.control['bc']['inlet']['p']
        self._pInBC = getattr(boundaryConditions, pBCDict['type'])(pBCDict)
        pBCDict = self.control['bc']['outlet']['p']
        self._pOutBC = getattr(boundaryConditions, pBCDict['type'])(pBCDict)
        # Loss factor BC
        zetaBCDict = self.control['bc']['inlet']['zeta']
        self._zetaInBC = getattr(boundaryConditions, zetaBCDict['type'])(zetaBCDict)
        zetaBCDict = self.control['bc']['outlet']['zeta']
        self._zetaOutBC = getattr(boundaryConditions, zetaBCDict['type'])(zetaBCDict)

        self.setInitialConditions(0)


    # Flow rate equation initial condition
    def setInitialConditions(self, t0):
        print('---> Starting Q0 iteration with initial Q0 = ', self.Q0)

        # Initialize the boundary conditions
        self._pOut = self._pOutBC.getValue(t0)
        self._pIn = self._pInBC.getValue(t0)
        self.Dp = self._pOut - self._pIn
        self._zetaOut = self._zetaOutBC.getValue(t0)
        self._zetaIn = self._zetaInBC.getValue(t0)
    
        # Initialize the flow rate
        L = -1
        convergence = [False, False]
        Qtol = 1E-10
        Q0old = self.Q0.copy()  # Old values, copied, otherwise we point to the same
        
        # for region in self.regions:
        #     region.update()

        self.constants(self.Q0)  # Initialize the flow constants
        for i in range(50):
            # print('Starting initial condition iteration: ', i + 1)
            for r, region in enumerate(self.regions):
                tio = 0.5 * self._fluid['rho'] * (self._zetaIn / region.data['s'][0] ** 2 +
                                                  self._zetaOut / region.data['s'][L] ** 2)
                tcv = self.Wc(r, 1, region)[L] + self.Wv(r, 1, region)[L]

                self.Q0[r] = np.sqrt(-self.Dp / (tio + tcv))  # New Q0

                self.constants(self.Q0)  # Constants for the new Q0

                if (np.abs(self.Q0[r] - Q0old[r])) / np.abs(self.Q0[r]) < Qtol:
                    convergence[r] = True

            Q0old = self.Q0.copy()

            # Check convergence of both regions
            if convergence[0] and convergence[1]:
                self.state = self.Q0
                print('---> Initial Q0 has converged to: ', self.Q0, '. Iterations: ', i)
                break

        if not convergence[0] and not convergence[1]:
            raise ValueError('     ERROR: The fluid initial condition has not converged!')
        
    # Update the region after a boundary modification
    def updateRegions(self):
        for i, region in enumerate(self.regions):
            region.update()
            
    def updatePressure(self, state):
        convergence = False
        ctol = 1E-3
        px = np.zeros(self.dof, dtype=object)
        Q = np.zeros(self.dof, dtype=object)
        dQ = np.zeros(self.dof, dtype=object)

        Q0 = state
        dQ0 = self.rhs
        
        for i, region in enumerate(self.regions):
            # Update the flow rate
            Q[i] = Q0[i] - region.data['dsi']
            
            # Update the acceleration
            dQ[i] = dQ0[i] - region.data['ddsi']

            # Correct the pressure for the new acceleration
            Wcv = (self.Wc(i, Q[i] ** 2, region) + self.Wv(i, Q[i] ** 2, region))
            px[i] = (self._pIn
                          - self._fluid['rho'] * 0.5 * self._zetaIn * (Q0[i] / region.data['s'][0]) ** 2
                          - self.Wt(i, dQ[i], region)
                          - Wcv)
            # Update the boundary pressure
            for boundary in region.boundaries():
                if boundary.isFlexible():
                    boundary.update("Forces", px[i])
            
            # Some debugging info
            #pnorm = np.linalg.norm(self.px[i] - px[i]) / np.linalg.norm(px[i])
            #print(" Pressure change norm is ", pnorm)
            #print("Pressure term:", self._pIn)
            #print("Convective-Viscous term:", Wcv)
            #print("Boundary term:", self._fluid['rho'] * 0.5 * self._zetaIn * (Q0[i] / region.data['s'][0]) ** 2)
            #print("Transient term:", self.Wt(i, dQ[i], region))
            
        return False
        

    # Update the state variable Q0, its derivative and the pressure after convergence.
    def advance(self, tf, state):
        
        # Update the flow rate and flow rate speed (last calculated RHS)
        self.Q0 = state
        self.state = state
        self.dQ0 = self.rhs
        
        for i, region in enumerate(self.regions):
            # Update the flow rate
            self.Q[i] = self.Q0[i] - region.data['dsi']
            self.vx[i] = self.Q[i] / region.data['s']

            # Update the acceleration
            self.dQ[i] = self.dQ0[i] - region.data['ddsi']

            # Correct the pressure for the new acceleration
            Wcv = (self.Wc(i, self.Q[i] ** 2, region) + self.Wv(i, self.Q[i] ** 2, region))
            self.px[i] = (self._pIn
                          - self._fluid['rho'] * 0.5 * self._zetaIn * (self.Q0[i] / region.data['s'][0]) ** 2
                          - self.Wt(i, self.dQ[i], region)
                          - Wcv)
            
            # Update the boundary pressure
            for boundary in region.boundaries():
                if boundary.isFlexible():
                    boundary.update("Forces", self.px[i])
            
            # Calculate the fluid power
            # Note that pi is not equal to p_In since the inlet is reservoir. This balance is calculated in terms of 
            # px[x] and not in terms of pIn (which would imply to take zero velocity)
            # pi, Qi = self.px[i][0], self.Q[i][0]
            # po, Qo = self.px[i][-1], self.Q[i][-1]
            # vi, vo = Qi / region.data['s'][0], Qo / region.data['s'][-1]
            
            # # This should be correct but it does not include the viscous dissipation
            # self.power = self.power.fromkeys(self.power, 0.0)
            # dpdx = np.gradient(self.px[i], self._mesh.x, edge_order=2)

            # self.power['inertial'] += self._fluid['rho'] * si.simps(self.dQ[i] * self.vx[i], self._mesh.x)
            # self.power['convective'] +=  0.5 * self._fluid['rho'] * (Qo * vo**2 - Qi * vi**2)
            # self.power['convective2'] +=  0.5 * self._fluid['rho'] * (
            #     dpdx[-1]**3 * region.data['s'][-1]**7 / (1120 * self._fluid['mu']**3) - 
            #     dpdx[0]**3 * region.data['s'][0]**7 / (1120 * self._fluid['mu']**3)) 
            # self.power['pressure'] +=  (Qo * po - Qi * pi)     
            # self.power['viscous'] += 0.0
            # fvisc = dpdx**2 * region.data['s']**5 / (120 * self._fluid['mu'])
            # self.power['heat'] += si.simps(fvisc, self._mesh.x)

        # Update the boundary condition
        self._pOut = self._pOutBC.getValue(tf)
        self._pIn = self._pInBC.getValue(tf)
        self.Dp = self._pOut - self._pIn
        self._zetaOut = self._zetaOutBC.getValue(tf)
        self._zetaIn = self._zetaInBC.getValue(tf)
        
    
    def initialize(self):
        # This method builds operators that are constant through the ODE solution (solid based operators)
        # This avoids the re-evaluation of the operators in the RHS calculation. The only downside is that 
        # the Fanning friction factor f is evaluated explicitly. 
        
        # Initialize the operators and temporary variables
        self._Wt_1_L = np.zeros(self.dof)
        self._t0 = np.zeros(self.dof)
        self._t1 = np.zeros(self.dof)
        self._t2 = np.zeros(self.dof)
        
        # Aliases
        rho = self._fluid['rho']
        zi = self._zetaIn
        zo = self._zetaOut

        
        L = -1
        for i, region in enumerate(self.regions):
            s = region.data['s']
            dsi = region.data['dsi']
            self._Wt_1_L[i] = self.Wt(i, 1, region)[L]
            Wt_ddsi_L = self.Wt(i, region.data['ddsi'], region)[L]
            Wcv_1_L = self.Wc(i, 1, region)[L] + self.Wv(i, 1, region)[L]
            temp = 2 * dsi
            Wcv_2dsi_L = self.Wc(i, temp, region)[L] + self.Wv(i, temp, region)[L]
            temp = dsi**2
            Wcv_dsi2_L = self.Wc(i, temp, region)[L] + self.Wv(i, temp, region)[L]
            
            self._t2[i] = 0.5 * rho * ( (zo / s[L]**2) + ( zi / s[0]**2) ) +  Wcv_1_L
            self._t1[i] = -rho * zo * dsi[L] / s[L]**2 - Wcv_2dsi_L
            self._t0[i] = 0.5 * rho * zo * dsi[L]**2 / s[L]**2 - Wt_ddsi_L + Wcv_dsi2_L
            
# Evaluate the RHS ofthe equation
    def RHS_old(self, time, state):
        
        # Update the state
        Q0 = state  # Update the flow rate
        Q = np.zeros(self.dof, dtype=object)
        self._ti = time

        L = -1  # Index of the end of the channel

        self.constants(Q0) # Calculate the constants
        
        # Assemble the rhs of each region
        self.rhs = np.zeros(self.dof)
        for i, region in enumerate(self.regions):            
            # Aliases
            s = region.data['s']
            dsi = region.data['dsi']
            
            # Flow rate function
            Q[i] = Q0[i] - dsi  
            
            # Calculate intermetdiate terms
            t2 = 0.5 * self._fluid['rho'] * (self._zetaIn * (Q[i][0] / s[0]) ** 2
                 + self._zetaOut * (Q[i][L] / s[L]) ** 2)
            Qi2 = Q[i]**2
            Wcv = (self.Wc(i, Qi2, region) + self.Wv(i, Qi2, region))
            t4 = - self.Wt(i, region.data['ddsi'], region)[L]  # Old value of acceleration
            
            self.rhs[i] = -(1 / self.Wt(i, 1, region)[L]) * (self.Dp + t2 + Wcv[L] + t4)
    
        return self.rhs


    # Evaluate the RHS ofthe equation
    def RHS(self, time, state):
        
        # Update the state
        Q0 = state 
        
        # Calculate the constants with the new state       
        self.constants(Q0) 

        # Assemble the rhs of each region
        self.rhs = np.zeros(self.dof)
        for i, region in enumerate(self.regions):            
            self.rhs[i] = - (self.Dp + self._t2[i] * Q0[i]**2 + self._t1[i] * Q0[i]  + self._t0[i]) / self._Wt_1_L[i]  

        return self.rhs


    # ----- Flow Operators ----- #
    # Transient Operator
    def Wt(self, i, fx, region):
        size = region.data['s']
        Wt = self._fluid['rho'] * si.cumtrapz(fx / size, self._mesh.x, initial=0)
        return Wt

    # Convective Operator
    def Wc(self, i, fx, region):
        size = region.data['s']
        wc = si.cumtrapz( (1/size) * np.gradient(fx / size, self._mesh.x, edge_order=2), self._mesh.x, initial=0)
        Wc = self._fluid['rho'] * self.xix[i] * wc
        return Wc

    # Viscous Operator
    def Wv(self, i, fx, region):
        size = region.data['s']
        Wv = 0.25 * self.f0[i] * si.cumtrapz(fx / size**3, self._mesh.x, initial=0)
        return Wv
    

    # Calculate some constants
    def constants(self, Q0):
        
        # Update the reference diameter and use the current flow rate before calculate the numbers
        for i, region in enumerate(self.regions):
            self.dRef[i] = region.data['s'][0]  # Reference channel size
            self.vRef[i] = Q0[i] / self.dRef[i]
            self.eRef[i] = self.dRef[i] / self.lRef
        
        self.calcNumbers()
        
        for i, region in enumerate(self.regions):
            self.Rd[i] = ReynoldsNumber(self, L=self.dRef[i], V=self.vRef[i])
            if self.Rd[i].value < 1000:  # Laminar
                if self.Rd[i].value <= 1: # Set the no flow values to zero to avoid indetermination
                    self.f0[i] = 1E-15
                    self.xix[i] = 0
                    # self.eta = -0
                else: # Laminar values
                    self.f0[i] = 48.0 / self.Rd[i].value
                    self.xix[i] = 6.0 / 5.0
                    # self.eta = -self.f0 / Q0
            else:  # Turbulent
                self.f0[i] = 0.26 * self.Rd[i].value ** -0.24
                self.xix[i] = 1.0
                # self.eta = -(0.0624 * self.Rd.value ** -0.24) / Q0

    def calcNumbers(self):
        super().calcNumbers()