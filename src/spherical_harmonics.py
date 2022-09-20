from inputs import *

# Rotate and scale 3D object to match point cloud
def rotate_scale(vertices, pc):
    # Rotate
    pc_mean = np.mean(pc, axis=0)
    pc = pc - pc_mean
    vertices = vertices - pc_mean
    pc_cov = np.cov(pc.T)
    vertices_cov = np.cov(vertices.T)
    U, S, V = np.linalg.svd(pc_cov.dot(vertices_cov.T))
    R = U.dot(V)
    vertices = vertices.dot(R.T)
    # Scale
    pc_scale = np.sqrt(np.sum(pc**2)/np.sum(vertices**2))
    vertices = vertices * pc_scale
    return vertices

def rotate_vector(v, theta, axis):
    v = np.asarray(v)
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    R = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                  [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                  [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    return np.dot(R, v)

def read(file):
    if file[-3:] == 'stl':
        return read_stl(file)
    elif file[-3:] == 'ply':
        return read_ply(file)
    elif file[-3:] == 'obj':
        return read_obj(file)
        

def read_stl(file):
    mesh_data = mesh.Mesh.from_file(file)
    num_faces = len(mesh_data.vectors)
    vertices = []
    for i in range(num_faces):
        for j in range(3):
            vertices.append(mesh_data.vectors[i][j])
    vertices = np.unique(np.asarray(vertices), axis=0)
    return vertices


def read_ply(file):
    pcd = o3d.io.read_point_cloud(file)
    vertices = np.unique(np.asarray(pcd.points), axis=0)
    return vertices


# A fucntion to read vertices from obj file
def read_obj(file):
    vertices = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('v '):
                vertices.append([float(x) for x in line.split()[1:4]])
    return np.asarray(vertices)


def sph2cart(theta, phi, r=1):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def cart2sph(x, y, z):
    phi = np.arctan2(y, x); neg_id = np.where(phi < 0);  phi[neg_id] += 2*np.pi
    theta = np.arctan2(z, (np.sqrt(x**2 + y**2))); theta = np.pi/2 - theta
    r = np.sqrt(x**2+y**2+z**2)
    return theta, phi, r


def plot_3D_pc(ax, x, y, z, maps=None, title_txt=''):
    fmax, fmin = maps.max(), maps.min()
    fcolors = (maps - fmin)/(fmax - fmin)
    _scatter = ax.scatter(x, y, z, marker=".", linewidth=0.05, facecolor=cm.seismic(fcolors))
    plt.axis('off')


def basis_functions(max_n, thetas, phis):
    Y_mat = np.zeros([thetas.shape[0], (max_n+1)**2], dtype=complex)
    for n in range(max_n + 1):
        for m in range(-n, n+1):
            if m >= 0:
                Y_mat[:, n**2 + n + m] = special.sph_harm(abs(m), n, phis, thetas)
            else:
                Y_mat[:, n**2 + n + m] = ((-1)**m) * special.sph_harm(abs(m), n, phis, thetas).conj()
    return Y_mat


def shape_descriptors(coefs, max_n):
    Dr = np.zeros((max_n, ))
    for n in range(max_n):
        for m in range(-n, n+1):
            Dr[n] += coefs[n**2 + n + m].real**2 + coefs[n**2 + n + m].imag**2
    # Dr /= Dr[0]
    return Dr


def math_basis_functions(max_n, thetas, phis):
    Y_mat = np.zeros([thetas.shape[0], thetas.shape[1], (max_n+1)**2], dtype=complex)
    for n in range(max_n + 1):
        for m in range(-n, n+1):
            if m >= 0:
                Y_mat[:, :, n**2 + n + m] = special.sph_harm(abs(m), n, phis, thetas)
            else:
                Y_mat[:, :, n**2 + n + m] = ((-1)**m) * special.sph_harm(abs(m), n, phis, thetas).conj()
    return Y_mat


def math_basis_reconstruction(coefs, rec_thetas, rec_phis, max_n):
    rec_Y_mat = math_basis_functions(max_n, rec_thetas, rec_phis)
    rec_rs = np.zeros(rec_thetas.shape, dtype=complex)
    for n in range(max_n + 1):
        for m in range(-n, n+1):
            rec_rs += coefs[n**2 + n + m] * rec_Y_mat[:, :, n**2 + n + m]
    return rec_rs.real


def plot_math_rec(fig, rec_max_n, coefs, sub, name):
    rec_thetas, rec_phis = np.linspace(0, np.pi, sub), np.linspace(0, 2*np.pi, sub)
    grid_thetas, grid_phis = np.meshgrid(rec_thetas, rec_phis)

    for i, max_n in enumerate(rec_max_n):
        rec_rs = math_basis_reconstruction(coefs, grid_thetas, grid_phis, max_n)
        rec_x, rec_y, rec_z = sph2cart(grid_thetas, grid_phis, rec_rs)
        fmax, fmin = rec_rs.max(), rec_rs.min()
        fcolors = (rec_rs - fmin)/(fmax - fmin)
        ax = fig.add_subplot(111, projection='3d')
        #plt.plot(projection='3d')
        surf = ax.plot_surface(rec_x, rec_y, rec_z, facecolors=cm.seismic(fcolors))
        #plt.title('{}'.format(max_n), y=-0.1)
        plt.axis('off')
    plt.show()
    filename = join(out_dir, 'rec_{}.png'.format(name))
    plt.savefig(filename)


if __name__ == '__main__':
    Drs = []
    for i, input_file in enumerate(input_files):
        file = os.path.join(input_dir, input_file)
        v = read(file)
        v = v - np.mean(v, axis=0)

        if input_file.endswith('.obj'):
            v /= 4.4

        thetas, phis, rs = cart2sph(v[:, 0], v[:, 1], v[:, 2])
        fig = plt.figure()
        fig.set_size_inches(2, 10)
        #fig = plt.figure(figsize=plt.figaspect(1.))
        ax = fig.add_subplot(111, projection='3d')
        plot_3D_pc(ax, v[:, 0], v[:, 1], v[:, 2], rs, names[i])

        Y_mat = basis_functions(max_n, thetas, phis)
        coefs = (np.linalg.inv(np.transpose(Y_mat).dot(Y_mat))
                .dot(np.transpose(Y_mat))).dot(rs)
        Drs.append(shape_descriptors(coefs, max_n))
        plot_math_rec(fig, rec_max_n=rec_max_n, coefs=coefs, sub=150, name=names[i])


    # np.savetxt('Drs.csv', np.asarray(Drs), delimiter=',')
    #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    #rc('text', usetex=True)
    #rc('font', size=24)
    rc('font', size=20)

    plt.figure()
    for i in range(len(names)):
        plt.loglog(np.arange(len(Drs[i])-1), (Drs[i][1::]), colors[i], label=names[i])
                #label=names[i], color=colors[i//5], 
                #alpha=1-(i%5)/(len(names)/len(colors)))



    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.69, 8.27)
    
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

    # plt.grid(which='minor')
    plt.ylabel(r'Shape descriptor $(D_r)$')
    plt.xlabel(r'Harmonic number ($r$)')
    # Set y-axis limit
    plt.ylim(10**-7, 10**-1)
    plt.legend(loc='upper right', ncol=2, frameon=1, fontsize=20)
    plt.savefig(join(out_dir,'shape_descriptors.pdf'), bbox_inches='tight')
    plt.show()
