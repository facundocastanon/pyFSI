from curses import A_BOLD
import numpy as np

class generalizedAlpha:
    def __init__(self, solid, time):
        self._solid = solid
        self._time = time
        
        self.beta = 0.25
        self.gamma = 0.50
        self.rho = solid.execution()["solver"]['integrators']['solid']['rho']
        
        # Old idea to avoid instabilities
        # self.rhoDamp = solid.execution()["solver"]['integrators']['solid'].get('rhoDampingTime', 0.0)
    
    def advance(self, tspan):
        # Damp oscillatory solution at the begginning of the simulation by reducing rho (too rude to work)
        # if tspan[-1] <= self.rhoDamp:
        #     rho = 0.1
        # else:
        #     rho = self.rho
        
        rho = self.rho
            
        state = np.zeros(self._solid.sof)  # State space DOFs
        solid, time = self._solid, self._time
        beta, gamma = self.beta, self.gamma
        # Parameters
        dt = tspan[1] - tspan[0]
        
        Uold = solid.a
        Vold = solid.da
        Aold = solid.dda
        Fold = solid.F
        Fnew = solid.forces()
        
        
        alpham = (2 * rho - 1) / (rho + 1)
        alphaf = rho / (rho + 1)
        gamma = 0.5 + alphaf - alpham
        beta = 0.25 * (gamma + 0.5)**2
        af = 1 - alphaf
        am = 1 - alpham
        
        gammap = gamma / (beta * dt)
        betap = am / (dt**2 * beta * af)
        
        c1 = am / (dt**2 * beta)
        c2 = gamma * af / (beta * dt)
        c4 = (am - 2 * beta) / (2 * beta)
        c5 = (gamma * af - beta) / (beta)
        c6 = (gamma - 2 * beta) * af / (2 * beta)
        c7 = c2 / af
        c8 = (gamma - beta) / beta
        c9 = (gamma - 2 * beta) / (2 * beta)
        c10 = c1 / am
        c11 = c10 * dt
        c12 = (1 - 2 * beta) / (2 * beta)
   
        # Effective stiffness
        Keff = c1 * solid.M + c2 * solid.C + af * solid.K
        
        # Effective force
        Faf = af * Fnew + alphaf * Fold
        FK = np.dot(solid.K, alphaf * Uold)
        FM = np.dot( solid.M, (c1 * Uold + c1 * dt * Vold + c4 * Aold))
        FC = np.dot(solid.C, (c2 * Uold + c5 * Vold + c6 * dt * Aold ))
        Feff = Faf - FK + FM + FC
        
        Unew = np.linalg.solve(Keff, Feff)
        
        DU = Unew - Uold
        Vnew = c7 * DU - c8 * Vold - c9 * dt * Aold 
        Anew = c10 * DU - c11 * Vold - c12 * Aold
        
        # Update the state
        state[0:solid.dof] = Unew
        state[solid.dof:solid.sof] = Vnew
        dda = Anew
        
        return state, dda