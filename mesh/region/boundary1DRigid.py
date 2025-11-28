from pyFSI.mesh.region.boundary1D import *

class boundary1DRigid(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DRigid'

    def __init__(self, execution, control, mesh):
        super().__init__(execution, control, mesh)
        
        self.dy = 0  # Obvious but neccessary for region calculations
        self.ddy = 0 # Obvious but neccessary for region calculations

        
    def isFlexible(self):  # The boundary is rigid
        return False

