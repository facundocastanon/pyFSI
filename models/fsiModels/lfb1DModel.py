# Change the system matrix.
# Corrected the position of G, Small acceleration and changed the sign of Tb
import copy, sys
import numpy as np, scipy.integrate as si
from pyFSI.post.prints import printArray

from pyFSI.vectors.eigen import eigenSystemVector as es
from pyFSI.models.fsiModels.fsiBase import fsiBase
from pyFSI.models.properties.dimensionlessNumbers import dimensionlessNumber


class lfb1D(fsiBase):
    def __repr__(self):
        return 'lfb1D'

    def __init__(self, execution, control, TIME, OREG):
        """Linear stability two-region leakage flow model

        Args:
            execution (dict): Execution control dictionary
            control (dict): Control dictionary
            TIME (object): Time object
            OREG (object): Object Registry
        """
        
        # Call the base class initialization
        super().__init__(execution, control, TIME, OREG)
        
        # Debug
        if self._debug:
            print("     WARNING: Region names are hardcoded.")

        # Size of the boundary state-space eigen matrix
        esize = self._flow._size # Size of the associated eigensystem
        self._esize = esize
        
        # Size of the FSI state-space system for two regions
        self._size = self._esize * 2 + 2

        # Initialize the system matrices
        self.K = np.zeros((esize, esize))  # Mass
        self.C = np.zeros((esize, esize))  # Damping
        self.M = np.zeros((esize, esize))  # Stiffness
        self.Tb = np.zeros(esize)
        self.Tt = np.zeros(esize)
        self.Bb = np.zeros(esize)
        self.Bt = np.zeros(esize)
        self.Db = np.zeros(esize)
        self.Dt = np.zeros(esize)
        self.Eb = np.zeros(esize)
        self.Et = np.zeros(esize)
        self.Gt = None
        self.Gb = None
        
        # Initialize the eigensystem for correct output.
        # Note that the eigenvector dimensions is not correct, but it is not important 
        # sice we re-create the object in solve)
        valuesComplex = np.zeros(esize, dtype=complex) + 1j
        vectorsComplex = np.zeros((esize, esize), dtype=complex) + 1j
        self.ES = es.eigenSystemVector.fromNumpy(valuesComplex, vectorsComplex)
        self.norm = np.zeros((esize, esize))  # Eigenvector Norm
        self.S = np.zeros((self._size, self._size))  # System matrix
        
        # Recontruction of the modal response
        self.tint = np.arange(0, 1E-1, 1E-3)
        self.shape = np.zeros((len(self.tint), self._solid.mesh().size))
        
        # Output variables
        self.varMap["eigenValues"] = "ES.values"
        self.varMap["eigenVectors"] = "ES.vectors"
        self.varMap["eigenValuesReal"] = "ES.valuesReal"
        self.varMap["eigenVectorsReal"] = "ES.vectorsReal"        
        self.varMap["eigenValuesImag"] = "ES.valuesImag"
        self.varMap["eigenVectorsImag"] = "ES.vectorsImag"
        self.varMap["shape"] = "shape"
        
    def solve(self, solver):
        # Note that it is not correct to calculate Y by solving eig(S.T) because
        # the eigenvectors are not normalized correctly. Normalization implies
        # Y.T = inv(X); this is the correct way to obtain a normalized system.
        
        # Solve the eigensystem
        self.R, self.X = np.linalg.eig(self.S)

        # Eigenvectors transpose
        self.Y = np.linalg.inv(self.X).T
        
        #  Create an eitenSystemVector object
        self.ES = es.eigenSystemVector.fromNumpy(self.R, self.X)
        
        # Compute the deformed shape
        if "shape" in self._execution['output'][self.name]:
            
            # Initial condition for selected modes of the state-space system.
            x0 = np.zeros(self._size)
            x0[self.control['responseModes']] = 1E-2
            
            # Calculate the evolution of the FSI generalized coordinates xt
            xt = solver.response(x0, self.tint, modes=self.control['responseModes'])
            
            # Reconstruct the transient motion of the solid as a linear combination 
            # of the FSI state-space evolution
            for i, t in enumerate(self.tint):
                self.shape[i, :] = self._solid.eigen.reconstruct(np.real(xt[0:self.ES.size, i]))

    def update(self):
        # Update the flow, the solid does not need to be updated
        self._flow.update()
        
        # Calculate dimensionless numbers
        self.calcNumbers()
        
        # Construct the FSI system matrix
        self.assemble()

    def assemble(self):
        # Aliases
        solid = self._solid
        flow = self._flow
        esize = self._esize

        # Calculate the norm
        # phiTphi = np.tensordot(solid.eigen.vectors, solid.eigen.vectors, axes=0)
        # # Construct the norm matrix (assume orthogonality of modes, only diagonal
        # # elements are calculated)
        # for i in range(esize):
        #     for j in range(esize):
        #         self.norm[i, j] = si.simps(phiTphi[i, j], solid.mesh().x)

        # Assemble the system matrices
        for i in range(0, esize):
            gi = solid.eigen.vectors[i]
            ki = flow.k[i]
            ci = solid.c[i] + flow.c[i]
            mi = solid.m[i] + flow.m[i]
            
            # Galerkin discretization
            for j in range(0, esize):
                gj = solid.eigen.vectors[j]
                self.K[j, i] = -si.simps(ki * gj, solid.mesh().x)
                self.C[j, i] = -si.simps(ci * gj, solid.mesh().x)
                self.M[j, i] = si.simps(mi * gj, solid.mesh().x)
        self.K -= solid.K
        
    
            # Assemble the system matrices
        # for i in range(0, esize):
        #     gi = solid.eigen.vectors[i]
        #     ki = solid.k[i] + flow.k[i]
        #     ci = solid.c[i] + flow.c[i]
        #     mi = solid.m[i] + flow.m[i]
        #     # Galerkin discretization
        #     for j in range(0, esize):
        #         gj = solid.eigen.vectors[j]
        #         self.K[j, i] = -si.simps(ki * gj, solid.mesh().x)
        #         self.C[j, i] = -si.simps(ci * gj, solid.mesh().x)
        #         self.M[j, i] = si.simps(mi * gj, solid.mesh().x)


        # System Region Vectors
        self.Gt = flow.Gq['channelTop']
        self.Gb = flow.Gq['channelBot']

        # Galerkin discretization
        # Note:
        for i in range(0, esize):
            gi = solid.eigen.vectors[i]
            self.Tt[i] = si.simps(flow.Tf['channelTop'] * gi, solid.mesh().x)  # Dim = n x 1
            self.Tb[i] = si.simps(flow.Tf['channelBot'] * gi, solid.mesh().x)  # Changed the sign # Dim = n x 1
            self.Bt[i] = si.simps(flow.Bq['channelTop'] * gi, solid.mesh().x)  # Dim = 1 x n
            self.Bb[i] = si.simps(flow.Bq['channelBot'] * gi, solid.mesh().x)  # Dim = 1 x n
            self.Dt[i] = si.simps(flow.Dq['channelTop'] * gi, solid.mesh().x)  # Dim = 1 x n
            self.Db[i] = si.simps(flow.Dq['channelBot'] * gi, solid.mesh().x)  # Dim = 1 x n
            self.Et[i] = si.simps(flow.Eq['channelTop'] * gi, solid.mesh().x)  # Dim = 1 x n
            self.Eb[i] = si.simps(flow.Eq['channelBot'] * gi, solid.mesh().x)  # Dim = 1 x n

        # Mass products
        Mi = np.linalg.inv(self.M)
        MiDotC = np.dot(Mi, self.C)
        MiDotK = np.dot(Mi, self.K)
        MiDotTb = np.dot(Mi, self.Tb)
        MiDotTt = np.dot(Mi, self.Tt)

        # System matrix
        self.S[0:esize, esize:2*esize] = np.identity(esize)
        self.S[esize:2*esize, 0:esize] = MiDotK
        self.S[esize:2*esize, esize:2*esize] = MiDotC
        self.S[esize:2*esize, 2*esize] = MiDotTb
        self.S[esize:2*esize, 2*esize+1] = -MiDotTt

        # Formulation considering acceleration terms
        if self.control['type'] == "Saravia":
            self.S[2 * esize, 0:esize] = -self.Eb - np.dot(self.Bb, MiDotK)
            self.S[2 * esize, esize:2 * esize] = -self.Db - np.dot(self.Bb, MiDotC)
            self.S[2 * esize, 2 * esize] = self.Gb - np.dot(self.Bb, MiDotTb)
            self.S[2 * esize, 2 * esize + 1] = np.dot(self.Bb, MiDotTt)

            self.S[2 * esize + 1, 0:esize] = self.Et + np.dot(self.Bt, MiDotK)
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt + np.dot(self.Bt, MiDotC)
            self.S[2 * esize + 1, 2 * esize] = np.dot(self.Bt, MiDotTb)
            self.S[2 * esize + 1, 2 * esize + 1] = self.Gt - np.dot(self.Bt, MiDotTt)

        # Formulation not considering acceleration terms
        elif self.control['type'] == "SaraviaReduced":
            self.S[2 * esize, 0:esize] = -self.Eb
            self.S[2 * esize, esize:2 * esize] = -self.Db
            self.S[2 * esize, 2 * esize] = self.Gb

            self.S[2 * esize + 1, 0:esize] = self.Et
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt
            self.S[2 * esize + 1, 2 * esize + 1] = self.Gt

        # Formulation of Tosi (I think some terms are wrong)
        elif self.control['type'] == "Tosi":
            self.S[2 * esize, 0:esize] = -(self.Eb + np.dot(self.Bb, np.dot(Mi, self.K)))
            self.S[2 * esize, esize:2 * esize] = -(self.Db + np.dot(self.Bb, np.dot(Mi, self.C)))
            self.S[2 * esize, 2 * esize] = -np.dot(self.Bb, np.dot(Mi, self.Tb))
            self.S[2 * esize, 2 * esize + 1] = self.Gb - np.dot(self.Bb, np.dot(Mi, self.Tt))

            self.S[2 * esize + 1, 0:esize] = self.Et + np.dot(self.Bt, np.dot(Mi, self.K))
            self.S[2 * esize + 1, esize:2 * esize] = self.Dt + np.dot(self.Bt, np.dot(Mi, self.C))
            self.S[2 * esize + 1, 2 * esize] = self.Gt + np.dot(self.Bt, np.dot(Mi, self.Tb))
            self.S[2 * esize + 1, 2 * esize + 1] = np.dot(self.Bt, np.dot(Mi, self.Tt))

        else:
            sys.exit("ERROR: No type in fsi formulation found...")
                
        #printArray(self.S, name="S")

    def calcNumbers(self):
        super().calcNumbers()
        self.dimNumbers['Mr'] = massRatio(self)
        self.dimNumbers['Kr'] = stiffnessRatio(self)
        self.dimNumbers['Gr'] = gapRatio(self)
        self.dimNumbers['Vp'] = viscousParameter(self)


# Dimensional numbers of this model
class massRatio(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Mr"
        mf = fsi.flow().fluid()['rho'] * fsi.flow().lRef  # Fluid mass
        ms = fsi.solid().material()['rho'] * fsi.solid().tRef  # Solid mass
        self.value = ms / mf


class stiffnessRatio(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Kr"
        kf = fsi.flow().fluid()['rho'] * fsi.flow().Q0**2 * fsi.flow().lRef**3 / fsi.flow().dRef**2
        ks = fsi.solid().material()['E'] * fsi.solid().control['I']   # Solid mass
        self.value = ks / kf


class gapRatio(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Gr"
        self.value = fsi.flow().eRef


class viscousParameter(dimensionlessNumber):
    def __init__(self, fsi):
        super().__init__()
        self.type = "Vp"
        self.value = fsi.flow().eRef**2 * fsi.flow().dimNumbers["Re"].value
