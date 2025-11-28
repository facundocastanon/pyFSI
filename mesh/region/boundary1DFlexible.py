from pyFSI.mesh.region.boundary1D import *

class boundary1DFlexible(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DFlexible'

    def __init__(self, execution, control, mesh):
        if "couple" in control:
            coupled = True
        else:
            coupled = False
        super().__init__(execution, control, mesh, coupled)

        # ----- Public attributes ----- #
        self.dy = np.zeros(mesh.size)   # Init to zero velocity
        self.ddy = np.zeros(mesh.size)  # Init to zero acceleration
        self.p = np.zeros(mesh.size)    # Boundary pressure function
        self.pc = 0                     # Pressure center
        self.f = np.zeros(mesh.size)    # Boundary Force

        # Output variable mapping
        self.varMap['pressure'] = 'p'
        self.varMap['pressureCenter'] = 'pc'
        self.varMap['totalForce'] = 'ft'

       
    def update(self, fieldName, fieldValues):
        if fieldName in self.control['couple']:           
            if fieldName == "Displacements":
                self.yOld = self.y
                self.y = self.y0 + fieldValues  # Set current position
            if fieldName == "Velocities":
                self.dyOld = self.dy
                self.dy = fieldValues  # Set current velocity
            if fieldName == "Accelerations":
                self.ddyOld = self.ddy
                self.ddy = fieldValues # Set current acceleration
            if fieldName == "Forces":
                self.pOld = self.p
                self.p = fieldValues
                self.fOld = self.f
                self.f = -self.p * self.n  # Current forces
                self.ft = si.simps(self.p, self.x)  
                self.pc = si.simps(self.p * self.x, self.x) / self.ft  # Pressure center

    def relax(self, fieldName, factor):
        if fieldName == "Displacements":
            self.y = self.yOld + (self.y - self.yOld) * factor  # Set current position
        if fieldName == "Velocities":
            self.dy = self.dyOld + (self.dy - self.dyOld) * factor  # Set current position
        if fieldName == "Accelerations":
            self.ddy = self.ddyOld + (self.ddy - self.ddyOld) * factor  # Set current position
        if fieldName == "Forces":
            self.f = self.fOld + (self.f - self.fOld) * factor  # Set current position

                    
    # ----- Abstract methods ----- #
    def isFlexible(self):
        return True
    
    def getPositions(self):
        return self.y
    
    def getDisplacements(self):
        return self.y - self.y0
    
    def getVelocities(self):
        return self.dy
    
    def getAccelerations(self):
        return self.ddy
    
    def getForces(self):
        return self.f
    
    def getPressures(self):
        return self.p

    def oldPositions(self):
        return self.yOld
    
    def oldDisplacements(self):
        return self.yOld - self.y0
    
    def oldVelocities(self):
        return self.dyOld
    
    def oldAccelerations(self):
        return self.ddyOld
    
    def oldForces(self):
        return self.fOld
    
    def oldPressures(self):
        return self.pOld