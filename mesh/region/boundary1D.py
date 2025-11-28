# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: Boundary 1D Base class
#    I    #     return:
# --------------------------------------------------------------------------- #
# Notes:
#   Interface for 1D boundaries
#
#
# --------------------------------------------------------------------------- #
from abc import ABCMeta, abstractmethod
import numpy as np
import scipy.integrate as si
from scipy.interpolate import interp1d
from pyFSI.execution.utilities import setParameter


class boundary1D(metaclass=ABCMeta):

    def __init__(self, execution, control, mesh, coupled=False):
        
        # ----- Public attributes ----- #
        self.name = control['name']
        self.size = mesh.size
        self.x = mesh.x
        self.y0 = None # Initial location
        self.y = None  # Current location
        self.n = control['normal']
        self.control = control

        # Output control
        self.output = []
        if execution['paths']['flowPath']:
            self.path = execution['paths']['flowPath']  # Associated path
        self.varMap = {}
        
        # self.d1 = None  # First derivative
        # self.ix = None  # Indefinite integral
        # self.iL = None  # Definite integral
        # self.dy = None  # First time derivative
        # self.ddy = None  # Second time derivative
        # self.dyi = None  # Indefinite spatial integral of the first time derivative
        # self.ddyi = None  # Indefinite spatial integral of the second time deriv
        
        # ----- Private attributes ----- #
        self._isCoupled = coupled  # if the field is coupled with precice or other
        self._mesh = mesh

        
        # Flag for parametric analysis
        if 'parameter' in control:
            self.parametric = True
        else:
            self.parametric = False
        
        # Procedures
        self.initialize()
        
        
    def initialize(self):
        # Execute the coordinate initialization
        self.y0 = getattr(self, self.control['method'])(self.control, self._mesh)

        # Transform the initial boundary position
        if 'transform' in self.control:
            self._transform()

        # Current Position of the boundary
        self.y = self.y0  
        
    # Default implementation of parametrize
    def parametrize(self, time):
        # Ramp the parameter using the proportional factor
        setParameter(self.control, time)
        
        self.__init__(self.control, self._mesh)


    # Transform the boundary 
    def _transform(self):
        if 'translate' in self.control['transform']:
            self.y0 += self.control['transform']['translate']


    # ----- Abstract methods ----- #
    # @abstractmethod
    # def update(self):
    #     pass

    # @abstractmethod
    # def isFlexible(self):
    #     pass
    
    # @abstractmethod
    # def getPositions(self):
    #     pass
    
    # @abstractmethod
    # def getDisplacements(self):
    #     pass
    
    # @abstractmethod
    # def getVelocities(self):
    #     pass
    
    # @abstractmethod
    # def getForces(self):
    #     pass
    
    # @abstractmethod
    # def getPressures(self):
    #     pass   
    
    
    @staticmethod
    def fromLineByTwoPoints(control, mesh):
        hi = control['hi']
        hf = control['hf']
        xvalues = mesh.x
        yvalues = np.zeros(len(xvalues))
        slope = (hf - hi) / (xvalues[-1] - xvalues[0])
        for i, x in enumerate(xvalues):
            yvalues[i] = slope * x + hi
        
        return yvalues
    
    @staticmethod
    def fromPoints(control, mesh):
        points = np.array(control['points'])
        if 'interpolation' in control:
            interpKind = control['interpolation']
        else:
            interKind = 'linear'
        interpolator = interp1d(points[0, :], points[1, :], kind=interpKind)
        
        return interpolator(mesh.x)


    # ----- Getters ----- #
    def __call__(self):
        return self.y

    def mesh(self):
        return self._mesh
    
    def isCoupled(self):
        return self._isCoupled

    def __repr__(self):
        return 'boundary1D'

    def write(self):
        pass