# This file works with Pymesh only

import sys, pymesh

class ICOSPHERE:
    def __init__(self):
        pass
    
    @staticmethod
    def construct(filename, refinment = 0, diameter = 1.0):
        sphere_mesh = pymesh.generate_icosphere(diameter,
                                                [0.,0.,0.],
                                                refinement_order= int(refinment))
        pymesh.save_mesh(filename, sphere_mesh);
        return sphere_mesh

if __name__ == "__main__":
    if len(sys.argv) == 3:
        filename  = str(sys.argv[0])
        refinment = sys.argv[1]
        diameter  = sys.argv[2]
    elif len(sys.argv) == 2:
        filename  = str(sys.argv[0])
        refinment = sys.argv[1]
        diameter  = 1.0
    elif len(sys.argv) == 1:
        filename  = str(sys.argv[0])
        refinment = 1
        diameter  = 1.0
    else:
        print("Error: the input variables are not defined correctly.")
    
    #try:
    mesh = ICOSPHERE.construct(filename, refinment, diameter)
    #except:
    #    print("Not enough/correct arguments have been passed.")