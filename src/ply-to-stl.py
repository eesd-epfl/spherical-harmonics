import pymeshlab
from inputs import *

def doc():
    pymeshlab.print_filter_list()
    pymeshlab.print_filter_parameter_list('compute_normals_for_point_sets')
    pymeshlab.print_filter_parameter_list('generate_surface_reconstruction_screened_poisson')


def generate_mesh(inp, out):
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(inp)
    ms.compute_normals_for_point_sets()
    ms.generate_surface_reconstruction_screened_poisson()
    ms.save_current_mesh(out)
    
    

if __name__ == '__main__':
    for input_file in input_files:
        inp = os.path.join(input_dir, input_file)
        out = inp[:-3] + 'stl'
        generate_mesh(inp, out)
        
    
