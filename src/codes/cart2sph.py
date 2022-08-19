# Mahmoud S. Shaqfa

# This function to convert fromt he cartesian to spherical coordiantes
# However, to save time we only get the azimuth and elevation angles as we always assume r = 1 when mapping to Spherical Harmonics
# Reference: http://mathworld.wolfram.com/SphericalCoordinates.html
# Those functions have been programmed assuming that:
#   Theta: the polar (elevation angle) measured from the positive Z-axis
#   Phi: the Azimuthal angle (longitudinal) measured fromt he positive X-axis
#   All the angles are in radians here

import numpy as np

def cart2sph(x, y, z, centroid):
    x-= centroid[0]; y-= centroid[1]; z-= centroid[2];
    y_x = y/x;
    y_x[y_x == -0.] = 0. # To account for arctan(-0.) = Pi
    azimuth = np.arctan(y_x)
    polar = np.arccos(z / np.sqrt(x**2. + y**2. + z**2.))
    return azimuth, polar

def sph2cart(azimuth, polar, r):
    x = r * np.sin(polar) * np.cos(azimuth)
    y = r * np.sin(polar) * np.sin(azimuth)
    z = r * np.cos(polar)
    return x, y, z

if __name__ == "__main__":
    # Test the function
    '''
    P3            P2
     _____________
    |             |      
    |             |
    |             |      
    |             |
    |             |      
    |_____________|

    P0           P1
    '''
    x = np.array([0., 1., 1., 0.])
    y = np.array([0., 0., 0., 0.])
    z = np.array([0., 0., 1., 1.])
    centroid =  np.array([0.5, 0.0, 0.5])
    az, el = cart2sph(x, y, z, centroid)
    print("The results| Azimuth: {},\n\t   Elevation: {}.".format(az, el))
    '''
         v3(3,0,4)
       *|
      * |
     *  |
    *   | v1(3,0,4)
      * |
        *  
         v2(4,0,2)
    '''
    x = np.array([2., 4., 3.])
    y = np.array([0., 0., 0.])
    z = np.array([3., 2., 4.])
    centroid =  np.array([0.0, 0.0, 0.0])
    az, el = cart2sph(x, y, z, centroid)
    print("The results| Azimuth: {},\n\t   Elevation: {}.".format(az, el))