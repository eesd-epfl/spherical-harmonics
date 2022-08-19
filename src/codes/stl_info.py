import numpy as np
from stl import mesh

class STLAnalyzer:
    def __init__(self, filename):
        self.filename = filename
        self.mesh_data = mesh.Mesh.from_file(self.filename)
        self.volume, self.cog, self.inertia = self.mesh_data.get_mass_properties()
        self.__analyze_surface()
    
    def __analyze_surface(self):
        self.num_faces = len(self.mesh_data.vectors)
        self.vertices = []
        for i in range(self.num_faces):
            # Three points to define a triangle (of a facet)
            for j in range(3):
                self.vertices.append(self.mesh_data.vectors[i][j])
        self.vertices = np.unique(np.asarray(self.vertices), axis=0) # clean duplicates
        self.num_verts = len(self.vertices)
        print("The number of faces is: {}, with {} vertices.".format(self.num_faces, self.num_verts))

if __name__ == "__main__":
    filename = 'test_stone.stl'
    test_stl = STLAnalyzer(filename)
