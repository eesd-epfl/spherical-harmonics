from stl_info import *

import numpy as np

class GeometryLoader:
    def __init__(self, filename):
        self.filename = filename
        self.mesh = STLAnalyzer(self.filename)
        self.vertices = self.mesh.vertices