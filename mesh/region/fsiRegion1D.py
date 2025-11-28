
import numpy as np
import scipy.integrate as si

class fsiRegion1D(object):

    def __repr__(self):
        return 'fsiRegion1D'

    def __init__(self, control, mesh, boundary):

        # ----- Public attributes ----- #
        self.name = control['name']
        self.type = control['type']

        # General container for derivatives, integrals, sizes
        self.data = {'s':       np.zeros(mesh.size),
                     'six':     np.zeros(mesh.size),  # Size indefinite integral
                     'siL':     0,  # Size definite integral between 0-L
                     'ds':      np.zeros(mesh.size),  # Vel of the size
                     'dds':     np.zeros(mesh.size),  # Vel of the size
                     'dsi':     np.zeros(mesh.size),  # Vel of the size integral
                     'ddsi':    np.zeros(mesh.size)}  # Accel of the size integral

        # ----- Private attributes ----- #
        self._mesh = mesh   # Reference to the mesh
        self._bBot = boundary[control['botBoundary']]  # Top boundary
        self._bTop = boundary[control['topBoundary']]   # Bottom boundary
        self._debug = mesh.debug()

        # ----- Procedures ----- #
        if self._debug:
            self.check()
            
        self.update()

    # Update the geometric data
    def update(self):
        x = self._mesh.x     
        
        self.data['s'] = self._bTop.y - self._bBot.y
        self.data['ds'] = self._bTop.dy - self._bBot.dy
        self.data['dds'] = self._bTop.ddy - self._bBot.ddy
        self.data['six'] = si.cumtrapz(self.data['s'], x, initial=0.0) 
        self.data['siL'] = si.simps(self.data['s'] , x) # Definite integral of the position
        self.data['dsi'] = si.cumtrapz(self.data['ds'], x, initial=0.0)
        self.data['ddsi'] = si.cumtrapz(self.data['dds'], x, initial=0.0)


    def check(self):
        # Check if the mesh density is ok
        if (self._mesh.x[1] - self._mesh.x[0]) > self.data['s'][0]:
            print("     WARNING: mesh dx is smaller the channel inlet "
                  "for region " + self.name + ". Results may be wrong. ")

    # ----- Getters ----- #
    # Return reference to the mesh
    def mesh(self):
        return self._mesh

    # Reference to the boundary objects
    def top(self):
        return self._bTop

    def bot(self):
        return self._bBot

    def boundaries(self):
        return [self._bTop, self._bBot]


