import numpy as np

class Newmark:
    def __init__(self, solid, time):
        self._solid = solid
        self._time = time
        
        self.beta = 0.25
        self.gamma = 0.50
    
    def advance(self, tspan):
        state = np.zeros(self._solid.sof)  # State space DOFs
        solid, time = self._solid, self._time
        beta, gamma = self.beta, self.gamma
        dt = tspan[1] - tspan[0]
        
        # Newmark parameters
        a0 = 1.0 / (beta * dt**2)
        a1 = gamma / (beta * dt)
        a2 = 1.0 / (beta * dt)
        a3 = 1.0 / (2.0 * beta) - 1.0
        a4 = gamma / beta - 1.0
        a5 = 0.5 * dt * (gamma / beta - 2.0)
        a6 = dt * (1.0 - gamma)
        a7 = gamma * dt  
        
        # Compute the residual
        KE = solid.K + a0 * solid.M + a1 * solid.C  # Effective Stiffness
        Rm = np.dot(solid.M, a0 * solid.a + a2 * solid.da + a3 * solid.dda)
        Rc = np.dot(solid.C, a1 * solid.a + a4 * solid.da + a5 * solid.dda)
        R = solid.forces() + Rm + Rc  # Effective force
        
        # print(" THE SOLID EFFECTIVE Rm NORM IS: ", np.linalg.norm(Rm))
        # print(" THE SOLID EFFECTIVE Rc NORM IS: ", np.linalg.norm(Rc))
        # print(" THE SOLID EFFECTIVE F NORM IS: ", np.linalg.norm(solid.forces()))
        # print(" THE SOLID EFFECTIVE FORCE NORM IS: ", np.linalg.norm(R))
        # Solve for the displacement, acceleration and velocity (respectively)
        
        u_new = np.linalg.solve(KE, R)
        ddu_new = a0 * (u_new - solid.a) - a2 * solid.da - a3 * solid.dda
        du_new = solid.da + a6 * solid.dda + a7 * ddu_new
        
        # Update the solid state (re-program, this is a very creapy design)
        state[0:self._solid.dof] = u_new
        state[self._solid.dof:self._solid.sof] = du_new
        dda = ddu_new
        
        return state, dda