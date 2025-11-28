import precice as prc
from colorama import Fore, Back, Style

from pyFSI.models.fsiModels.fsiBase import fsiBase

class preCICE(fsiBase):
    def __repr__(self):
        return 'precice_FSI'

    def __init__(self, execution, control, TIME, OREG):
        super().__init__(execution, control, TIME, OREG)
        
        # Select which is the object to solve
        # if control['type'] == 'SOLID':
        #     self.obj = solid
        # elif control['type'] == 'FLOW':
        #     self.obj = flow
        # elif control['type'] == 'MAGNETISM'
        #     self.obj = magnetism
        self.participant = OREG.get(self.control['participants'][0])


        # Precice initialization
        print("---> Initializing the precice interfaces...")
        self.interface = prc.Interface(self.participant.name, "../precice-config.xml", 0, 1)

        # Create the precice mesh objects (which are just pyFSI boundary objects)
        # Note that the flow has the mesh
        print("---> Creating the pyFSI-precice mesh objects...")
        self.meshes = []
        for mesh in self.control['meshes']:
            self.meshes.append(preciceMesh(mesh['name'],
                               mesh['read'],
                               mesh['write'],
                               self.interface,
                               OREG.get(mesh['boundary'])))

# A class to handle the mesh object in preCICE, actually it is a cloud of points
# on which the coupled variables are projected.
class preciceMesh:
    def __init__(self, name, read, write, interface, boundary):
        # Public Attributes
        self.name = name
        self.ID = interface.get_mesh_id(name)
        self.read = {}
        self.write = {}
        # Private attributes
        self._boundary = boundary # The preCICE mesh has a pyFSI boundary associated
        
        # Build the dictionary of IDs for the fields
        for field in read:
            self.read[field] = interface.get_data_id(field, self.ID)
        for field in write:
            self.write[field] = interface.get_data_id(field, self.ID)

        # Set the 3D fluid mesh coordinates (the 3D boundary coordinates)
        self.vertexIDs = interface.set_mesh_vertices(self.ID, boundary.positions)
        
        print('---> The preCICE mesh object named ', self.name, ' has been created succesfully.' )

    def readFields(self, interface):
        for fieldName, fieldID in self.read.items():
            fieldValues = interface.read_block_vector_data(fieldID, self.vertexIDs)
            print(Fore.BLUE + "----> Reading field " + fieldName + " for mesh " +  self.name)
            self._boundary.set3DFields(fieldName, fieldValues)  # Store the field in the boundary object

    def writeFields(self, interface, obj):
        for fieldName, fieldID in self.write.items():
            if obj.base == "flow":
                if fieldName == "Forces":
                    values = self._boundary.forces  # Pressure integral
                    print(Fore.BLUE + "----> Force written to the solid mesh")
            if obj.base == "solid":
                if fieldName == "Displacements":
                    values = self._boundary.displacements  # Pressure integral
                    print(Fore.BLUE + "----> Displacements written to the flow mesh")
                if fieldName == "Velocities":
                    
                    values = self._boundary.velocities  # Pressure integral
                    print(Fore.BLUE + "----> Velocities written to the flow mesh")
                if fieldName == "Accelerations":
                    values = self._boundary.accelerations  # Pressure integral
                    print(Fore.BLUE + "----> Accelerations written to the flow mesh")
            interface.write_block_vector_data(fieldID, self.vertexIDs, values)
    
    def boundary(self):
        return self._boundary