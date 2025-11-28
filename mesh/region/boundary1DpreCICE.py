from pyFSI.mesh.region.boundary1D import *

class boundary1DpreCICE(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DFlexible'

    def __init__(self, execution, control, mesh):
        super().__init__(execution, control, mesh, coupled=True)

        # ----- Public attributes ----- #
        self.dy = np.zeros(mesh.size)  # Init to zero velocity
        self.ddy = np.zeros(mesh.size) # Init to zero acceleration
        self.p = np.zeros(mesh.size) # Boundary Pressure
        self.f = np.zeros(mesh.size) # Boundary Force
        
        # Create the 3D Mesh for coupling
        self._make3DMesh(mesh)

          
    def update(self, fieldName, fieldValues):
        if fieldName == "Displacements":
            self.y = self.y0 + fieldValues  # Set current position
            self.positions[0:self.size, 1] = self.y
            self.positions[self.size:self.size*2, 1] = self.y
            self.displacements[0:self.size, 1] = fieldValues
            self.displacements[self.size:self.size*2, 1] = fieldValues
        if fieldName == "Velocities":
            self.dy = fieldValues  # Set current position
            self.velocities[0:self.size, 1] = self.dy
            self.velocities[self.size:self.size*2, 1] = self.dy   
        if fieldName == "Accelerations":
            self.ddy = fieldValues # Set current position
            self.accelerations[0:self.size, 1] = self.ddy
            self.accelerations[self.size:self.size*2, 1] = self.ddy           
        if fieldName == "Forces":
            self.p = fieldValues# Set current position
            self.f = -self.p * self.n
            self.forces[0:self.size, 1] = self.f
            self.forces[self.size:self.size*2, 1] = self.f
            
            
    # Update the 3D and the 1D fields using the preCICE adapter scheme
    def set3DFields(self, fieldName, fieldValues):
        if fieldName == "Displacements":
            self.y = self.y0 + fieldValues[0:self.size, 1]  # Set current position
            self.displacements = fieldValues
            self.positions = self.reference + self.displacements
            # print(fieldName)
            # print(fieldValues)
        if fieldName == "Velocities":
            self.dy = fieldValues[0:self.size, 1]  # Set current position
            self.velocities = fieldValues
        if fieldName == "Accelerations":
            self.ddy = fieldValues[0:self.size, 1]  # Set current position
            self.accelerations = fieldValues
        if fieldName == "Forces":
            #self.p = fieldValues[0:self.size, 1]  # Set current position
            self.f = fieldValues[0:self.size, 1]
            self.forces = fieldValues
            # print(fieldName)
            # print(fieldValues)
            

          
            
    def _make3DMesh(self, mesh):
        # Generate a cloud of points representing a 3D surface from the 1D boundary coordinates
        self.positions = np.zeros((mesh.size * 2, 3))  # The 3D mesh
        self.reference = np.zeros_like(self.positions) # Reference position
        self.displacements = np.zeros_like(self.positions)
        self.velocities = np.zeros_like(self.positions)
        self.accelerations = np.zeros_like(self.positions)
        self.forces = np.zeros_like(self.positions)
        
        # Set the initial positions
        self.positions[0:mesh.size, 0] = mesh.x
        self.positions[0:mesh.size, 1] = self.y
        self.positions[mesh.size:mesh.size*2, 0] = mesh.x
        self.positions[mesh.size:mesh.size*2, 1] = self.y
        self.positions[mesh.size:mesh.size*2, 2] = np.ones(mesh.size)  # The z=1 nodes        
        # self.Forces = np.zeros_like(self.vertices)  # Force to transfer

        
    # ----- Abstract methods ----- #
    def isFlexible(self):
        return True

