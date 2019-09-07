#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy
numpy.set_printoptions(precision=2, threshold=numpy.inf)

grid_size = 0.5

mesh_vertices = []
mesh_faces = []
def parseObj():
    for line in open("bread1.obj", 'r'):
        line_broken = line.rstrip().split(" ")
        if line_broken[0] == 'v':
            mesh_vertices.append([float(line_broken[1]), float(line_broken[2]), float(line_broken[3])])
        elif line_broken[0] == 'f':
            v0_index = int(line_broken[1].split("//")[0])-1
            v1_index = int(line_broken[2].split("//")[0])-1
            v2_index = int(line_broken[3].split("//")[0])-1
            mesh_faces.append([v0_index, v1_index, v2_index])

                
parseObj()
mesh_faces = numpy.array(mesh_faces, dtype=int)
mesh_vertices = numpy.array(mesh_vertices) # vertex positions
display_norms = numpy.zeros([len(mesh_faces) * 3, 3]) # face normals
display_verts = numpy.zeros([len(mesh_faces) * 3, 3]) # face verts
display_colrs = numpy.zeros([len(mesh_faces) * 3, 3]) # face verts colors

N_v = len(mesh_vertices)

verts_loc_grid = numpy.zeros([N_v, 3], dtype=int) # which grid the vertex is in
verts_loc_tet = numpy.zeros(N_v, dtype=int) # which tet the vertex is in
verts_color = numpy.zeros([N_v, 4]) # colors of vertices

mesh_itrpltn = numpy.zeros([N_v, 4]) # vertex position as interpolation of node values
mesh_min = numpy.zeros(3) # min of vertices in each dimension
mesh_max = numpy.zeros(3) # max of vertices in each dimension

grid_min = numpy.zeros(3, dtype=int)
grid_max = numpy.zeros(3, dtype=int)
grid = [] # whether each grid is occupied
grid2node = [] # grid (point) to node index
grid2tets = [] # grid to index of first tet

nodes = [] # positions of nodes

left_half_nodes = numpy.zeros(0) # nodes in the left half:

tetrahedra = [] # index of 4 nodes for each tet
bndry_nodes = [] # index of nodes on boundary


for i in range(0, 3):
    mesh_min[i] = numpy.amin(mesh_vertices[:,i])
    mesh_max[i] = numpy.amax(mesh_vertices[:,i])

grid_min = numpy.floor(mesh_min / grid_size).astype(int)
grid_max = numpy.ceil(mesh_max / grid_size).astype(int)
grid = numpy.zeros(grid_max - grid_min, dtype=int)
grid2node = numpy.zeros(grid_max + 1 - grid_min, dtype=int) - 1
grid2tets = numpy.zeros(grid_max - grid_min, dtype=int)


def createGrids():
    global nodes, grid, grid2node
    grid_cnt = 0
    # fill grid according to vertices positions
    for i in range (0, len(mesh_vertices)):
        vert = mesh_vertices[i]
        c0, c1, c2 = numpy.array( numpy.floor(vert / grid_size) - grid_min ).astype(int)
        c0, c1, c2 = numpy.minimum(numpy.array(grid.shape)-1, numpy.array([c0, c1, c2]))
        # cell does not exist, create nodes
        if not grid[c0, c1, c2] == 1:
            grid[c0, c1, c2] = 1
        verts_loc_grid[i] = numpy.array([c0, c1, c2])
    # fix the inside of the grid
    for c1 in range (0, grid.shape[1]):
        for c2 in range (0, grid.shape[2]):
            start_index = end_index = -1
            # find start_index
            for c0 in range (0, grid.shape[0]):
                if grid[c0, c1, c2] == 1:
                    start_index = c0
                    break
            # find end_index
            for c0 in reversed(range (0, grid.shape[0])):
                if grid[c0, c1, c2] == 1:
                    end_index = c0
                    break
            if start_index >= 0:
                for c0 in range (start_index, end_index+1):
                    if not grid[c0, c1, c2] == 1:
                        grid[c0, c1, c2] = 1
                    grid_cnt += 1
    print("grid_cnt: ", grid_cnt)
                        
def createNodes():
    global grid2node, left_half_nodes, nodes
    cnt = 0
    for c1 in range (0, grid.shape[1]):
        for c2 in range (0, grid.shape[2]):
            for c0 in range (0, grid.shape[0]):
                if grid[c0, c1, c2] > 0:
                    for i in [0, 1]:
                        for j in [0, 1]:
                            for k in [0, 1]:
                                n0, n1, n2 = numpy.array([c0, c1, c2]) + numpy.array([i, j, k])
                                # node does not exist, create node
                                if grid2node[n0, n1, n2] < 0:
                                    grid2node[n0, n1, n2] = cnt
                                    node_pos = grid_size * (numpy.array([n0, n1, n2]) + grid_min)
                                    nodes.append(node_pos)
                                    if node_pos[0] <= 0:
                                        left_half_nodes = numpy.append(left_half_nodes, cnt)                                   
                                    cnt += 1


def getBndryNodes():
    global bndry_nodes
    for j in range (0, grid.shape[1]):
        for k in range (0, grid.shape[2]):
            start_index = end_index = -1
            # find start_index
            for i in range (0, grid.shape[0]):
                if grid[i, j, k] > 0:
                    start_index = i
                    break
            # find end_index
            for i in reversed(range (0, grid.shape[0])):
                if grid[i, j, k] > 0:
                    end_index = i
                    break
            if start_index >= 0:
                for m in [0, 1]:
                    for n in [0, 1]:
                        bndry_nodes.append(grid2node[start_index, j+m, k+n])
                        bndry_nodes.append(grid2node[end_index+1, j+m, k+n])
    for i in range (0, grid.shape[0]):
        for k in range (0, grid.shape[2]):
            start_index = end_index = -1
            # find start_index
            for j in range (0, grid.shape[1]):
                if grid[i, j, k] > 0:
                    start_index = j
                    break
            # find end_index
            for j in reversed(range (0, grid.shape[1])):
                if grid[i, j, k] > 0:
                    end_index = j
                    break
            if start_index >= 0:
                for m in [0, 1]:
                    for n in [0, 1]:
                        bndry_nodes.append(grid2node[i+m, start_index, k+n])
                        bndry_nodes.append(grid2node[i+m, end_index+1, k+n])
    for i in range (0, grid.shape[0]):
        for j in range (0, grid.shape[1]):
            start_index = end_index = -1
            # find start_index
            for k in range (0, grid.shape[2]):
                if grid[i, j, k] > 0:
                    start_index = k
                    break
            # find end_index
            for k in reversed(range (0, grid.shape[2])):
                if grid[i, j, k] > 0:
                    end_index = k
                    break
            if start_index >= 0:
                for m in [0, 1]:
                    for n in [0, 1]:
                        bndry_nodes.append(grid2node[i+m, j+n, start_index])
                        bndry_nodes.append(grid2node[i+m, j+n, end_index+1])
    bndry_nodes = numpy.unique(bndry_nodes) 
    

def createTets():
    cnt = 0
    global grid, grid2node, tetrahedra
    for i in range (0, grid.shape[0]):
        for j in range (0, grid.shape[1]):
            for k in range (0, grid.shape[2]):
                if grid[i, j, k] == 1:
                    origin = grid2node[i][j][k]
                    back = grid2node[i][j][k+1]
                    right = grid2node[i+1][j][k]
                    back_right = grid2node[i+1][j][k+1]
                    up = grid2node[i][j+1][k]
                    up_back = grid2node[i][j+1][k+1]
                    up_right = grid2node[i+1][j+1][k]
                    up_right_back = grid2node[i+1][j+1][k+1]

                    grid2tets[i, j, k] = cnt
                    for m in range (0, 5):
                        tetrahedra.append(numpy.zeros([4]))
                    
                    tetrahedra[cnt] = numpy.array([origin, back_right, back, up_back])
                    tetrahedra[cnt + 1] = numpy.array([right, back_right, origin, up_right])
                    tetrahedra[cnt + 2] = numpy.array([up_back, up_right, up, origin])
                    tetrahedra[cnt + 3] = numpy.array([up_right_back, up_right, up_back, back_right])
                    tetrahedra[cnt + 4] = numpy.array([up_back, up_right, back_right, origin])
                    cnt += 5
    tetrahedra = numpy.array(tetrahedra)
    print("tetrahedra: ", tetrahedra.shape)
                    
def getIntrpltn():
    for i in range (0, len(mesh_vertices)):
        vert = mesh_vertices[i]
        n1, n2, n3 = verts_loc_grid[i]
        tet0_index = grid2tets[n1, n2, n3]
        for j in range (0, 5):
            tet = tetrahedra[tet0_index + j]
            # calculate barycentric coordinate
            col1 = nodes[tet[0]] - nodes[tet[3]]
            col2 = nodes[tet[1]] - nodes[tet[3]]
            col3 = nodes[tet[2]] - nodes[tet[3]]
            T = (numpy.matrix([col1, col2, col3])).T
            coord = numpy.array( numpy.matmul(numpy.linalg.inv(T), (vert - nodes[tet[3]])) )
            coord = numpy.array(coord[0])
            coord = numpy.append(coord, 1-numpy.sum(coord))
            if coord[0] >= 0 and coord[1] >= 0 and coord[2] >= 0:
                mesh_itrpltn[i] = coord
                verts_loc_tet[i] = tet0_index + j
                break
                
            
createGrids()
createNodes()
nodes = numpy.array(nodes)
print("nodes: ", nodes.shape)
createTets()                  
tetrahedra = numpy.array(tetrahedra)
getIntrpltn()

getBndryNodes()
# print("bndry_nodes: ", bndry_nodes)


# # The physics

# In[2]:


# number of nodes
N_n = nodes.shape[0]

# number of tetrahedra
N_t = tetrahedra.shape[0]

# damping constant
damp_constant = 5

# the displacements
u = numpy.zeros([N_n, 3])

# the velocities
v = numpy.zeros([N_n, 3])

# the force on each node
F = numpy.zeros([N_n, 3])

# the force resulting from heat
h = numpy.zeros([N_n, 3])

# B matrices for stress
B_s = numpy.zeros([N_t, 6, 12])

# B_s transpose multiplied by D
Bs_T_D = numpy.zeros([N_t, 12, 6])

# entries of shape function matrix b, c, d
coeffs = numpy.zeros([N_t, 3, 4])

# stiffness matrix
K_s = numpy.zeros([int(N_n*3), int(N_n*3)])
 
# gravity
g = numpy.array([0, -9.8, 0])

# mass of each node
m = 1 #kg

# external force
F_ex = numpy.zeros([N_n, 3])

# modulus of elasticity
E = 1e2

# poisson ratio
p_r = 0.3

# floor y position
floor_y = -1

# rotation angle
ry = 0

# time step, don't use actual dt which is not stable
dt = 0.01

# initialize D-matrix for Hooke's Law
D = numpy.identity(6)
D[0][0] = D[1][1] = D[2][2] = 1 - p_r
D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = p_r
D[3][3] = D[4][4] = D[5][5] = (1 - 2 * p_r) / 2
D *= E/((1 + p_r) * (1 - 2 * p_r))

# calculate N matrix, get its entries
def getCoeffs():
    global coeffs, tetrahedra, nodes
    for j in range (0, N_t): 
        t = tetrahedra[j]
        A = numpy.identity(4)
        for i in range (0, 4):
            A[i] = numpy.array([1,nodes[t[i]][0],nodes[t[i]][1],nodes[t[i]][2]])
        J = numpy.linalg.det(A)
        for k in range (1, 4):
            M = numpy.delete(A, k, axis=1)
            for i in range (0, 4):
                coeffs[j][k-1][i] = numpy.linalg.det(numpy.delete(M, i, axis=0))
                if (k + i) % 2 == 1:
                    coeffs[j][k-1][i] *= (-1)
        coeffs[j] /= J
        
# calculate B matrix for tets
def updateBsMatrices():
    global B_s, coeffs
    for j in range (0, N_t):
        for k in range (0, 3):
            for i in range (0, 4):
                B_s[j][k][i * 3 + k] = coeffs[j][k][i]
        for i in range (0, 4):
            B_s[j][3][i * 3] = 0.5 * coeffs[j][1][i]
            B_s[j][3][i * 3 + 1] = 0.5 * coeffs[j][0][i]
        for i in range (0, 4):
            B_s[j][4][i * 3] = 0.5 * coeffs[j][2][i]
            B_s[j][4][i * 3 + 2] = 0.5 * coeffs[j][0][i]
        for i in range (0, 4):
            B_s[j][5][i * 3 + 1] = 0.5 * coeffs[j][2][i]
            B_s[j][5][i * 3 + 2] = 0.5 * coeffs[j][1][i]  

# assemble stiffness matrix
def updateKs():
    global K_s, B_s, Bs_T_D
    for i in range (0, N_t):
        Bs_T_D[i] = numpy.matmul(B_s[i].T, D)
        K_l = numpy.matmul(Bs_T_D[i], B_s[i])
        S_1 = numpy.zeros([12, N_n * 3])
        for j in range (0, 4):
            for k in range (0, 3):
                S_1[:, 3 * tetrahedra[i][j] + k] = K_l[:,3 * j + k]
        for j in range (0, 4):
            for k in range (0, 3):
                K_s[3 * tetrahedra[i][j] + k] += S_1[3 * j + k]

# apply impact
def colDtctn():
    global F_ex, v, u
    F_ex.fill(0) # reset external force
    F_ex += g
    for i in range(0, N_n):
        node = nodes[i]
        if node[1] + u[i][1] <= floor_y:
#             print("reach floor")
            u[i][1] = floor_y - node[1]
            v[i][1] = - 0.5 * v[i][1] # apply some damping

def updateForce(dt_actual):
    global u, v, F, K, dt, F_ex, h
    F.fill(0) # reset force
    colDtctn() # do collision detection to get degrees of freedom
    F_in = numpy.array( (numpy.matmul(K_s, u.flatten())).reshape([N_n,3])) # calculate internal force
    F += F_ex # add external force
    F -= v * damp_constant # account for damping
    F += h # add heat force
    F -= F_in # minus internal force
    v += dt * F / m
#     print("v: ", v)    
    u += dt * v
#     print ("u: ", u)

getCoeffs()  
updateBsMatrices()
updateKs()

# def moveDown(dt_actual):
#     global u
#     u += numpy.array([0, -0.1, 0])
# moveDown(0)


# In[ ]:


# conductivity of material
k_cndct = 2.36e-1

# strain-temperature constant
alpha = 1e-2

# heat capacitance of material
c = 0.9

# heat radiation constant
ems_sig = 2.5e-1

# initial temperature
T_init = 0 #celcius

# env temperature
T_env = 100 # celcius

# each node's temperature
T = T_init * numpy.ones(N_n)

# force due to heat
h = numpy.zeros([N_n, 3])

# each node's heat flux
q_ex = numpy.zeros(N_n)

# B matrix for heat for each tet
B_h = numpy.zeros([N_t, 3, 4])

# stiffness matrix for heat
K_h = numpy.zeros([N_n, N_n])

# calculate B matrix for heat for tets
def updateBhMatrices():
    global B_h
    for j in range (0, N_t): 
        B_h[j][0] = coeffs[j][0]

# assemble stiffness matrix for heat
def updateKh():
    global K_h, B_h, k_cndct
    for i in range (0, N_t):
        K_e = k_cndct * numpy.matmul(B_h[i].T, B_h[i])
        for j in range (0, 4):
            for k in range (0, 4):
                K_h[tetrahedra[i][j]][tetrahedra[i][k]] += K_e[j][k]
    
# calculate heat exchange on boundary nodes              
def getBndryFlux():
    global T, bndry_nodes
    bndry_flux = numpy.zeros(N_n)
    for i in bndry_nodes:
        bndry_flux[i] = ems_sig * ( T_env - T[i] )
    return bndry_flux


# calculate "heat" force
def updateh():
    global Bs_T_D, T, h
    h.fill(0)
    for i in range (0, N_t):
        avg_T = ( T[tetrahedra[i][0]] + T[tetrahedra[i][1]] + T[tetrahedra[i][2]] + T[tetrahedra[i][3]] ) / 4.0
        strain_h = numpy.array([alpha * avg_T, alpha * avg_T, alpha * avg_T, 0, 0, 0])
        h_l = numpy.matmul(Bs_T_D[i], strain_h)
        for j in range (0, 4):
            h[tetrahedra[i][j]][0] += h_l[j*3]
            h[tetrahedra[i][j]][1] += h_l[j*3+1]
            h[tetrahedra[i][j]][2] += h_l[j*3+2]

def updateHeat(dt_actual):
    global c, m, T, K_h, dt
    q_in = - numpy.matmul(K_h, T) # internal heat transfer
    q_ex = getBndryFlux() # external heat transfer
    q = q_in + q_ex
    dQ = q * dt
    T += (1.0 / (c * m)) * dQ 
    
def animate(dt_actual):
    global ry
    updateHeat(dt_actual)
    updateh()
    updateForce(dt_actual)
    ry += 0.2
    
updateBhMatrices()
updateKh()


# # render 

# In[ ]:


from pyglet.gl import *
import pyglet

gl_vertices = list()
gl_normals = list()
gl_colors = list()

# number of vertices
verts_num = 0

# number of faces
N_f = 0

# render surface only
render_surf = False

# render half of cube
render_half = False

# save images
save_imgs = False
frame_num = 1


window = pyglet.window.Window()

def setup():
    glClearColor(0, 0, 0, 0)
    glEnable(GL_DEPTH_TEST)
    if not render_surf:
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glEnable(GL_LIGHT1)

    # Define a simple function to create ctypes arrays of floats:
    def vec(*args):
        return (GLfloat * len(args))(*args)
    
    glLightfv(GL_LIGHT0, GL_POSITION, vec(20, 30, 20))
#     glLightfv(GL_LIGHT0, GL_SPOT_EXPONENT, (GLfloat)(1))
#     glLightfv(GL_LIGHT0, GL_SPECULAR, vec(0.2, 0.2, 0.2, 1))
    glLightfv(GL_LIGHT0, GL_DIFFUSE, vec(1.0, 1.0, 1.0, 1))
    
#     glLightfv(GL_LIGHT1, GL_POSITION, vec(0, 0, 5))
#     glLightfv(GL_LIGHT1, GL_DIFFUSE, vec(.5, .5, .5, 1))
#     glLightfv(GL_LIGHT1, GL_SPECULAR, vec(1, 1, 1, 1))
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, vec(193./255., 217./255., 1.0, 0.5))
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, vec(1, 1, 1, 0.5))
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50)

def prepareVerts():
    global gl_vertices, gl_normals, gl_colors, verts_num, N_f
    indices = [[1, 2, 0], 
             [1, 3, 2], 
             [2, 3, 0], 
             [0, 3, 1]]
    normals = []
    vertices = []
    colors = []
    nodes_pos = nodes + u
    def inLeftHalf(j, i):
        for k in range (0,3):
            if not (tetrahedra[j][indices[i][k]] in left_half_nodes):
                return False
        return True
    for j in range (0, tetrahedra.shape[0]):
        for i in range (0, 4):
            if render_half and not inLeftHalf(j, i):
                continue
            n = numpy.cross(nodes_pos[tetrahedra[j][indices[i][2]]] - nodes_pos[tetrahedra[j][indices[i][1]]], 
                            nodes_pos[tetrahedra[j][indices[i][0]]] - nodes_pos[tetrahedra[j][indices[i][1]]])
            n = n / numpy.linalg.norm(n)
            normals.extend(n.tolist() * 3)
            for k in range (0, 3):
                vertices.extend(nodes_pos[tetrahedra[j][indices[i][k]]])
#                 colors.extend(to_tmp_color(T[tetrahedra[j][indices[i][k]]]))
#                 colors.extend(to_bndry_color(tetrahedra[j][indices[i][k]]))
                colors.extend(0.8 * numpy.ones(4))
    verts_num = len(vertices)
    N_f = int(verts_num / 9)
    gl_vertices = (GLfloat * len(vertices))(*vertices)
    gl_normals = (GLfloat * len(normals))(*normals) 
    gl_colors = (GLfloat * len(colors))(*colors) 

def getMeshVertices():
    nodes_pos = nodes + u
    global verts_loc_tet, mesh_itrpltn, mesh_vertices, N_f
#     mesh_vertices.fill(0)
    verts_color.fill(0)
    for i in range (0, N_v):
        tet = tetrahedra[verts_loc_tet[i]]
        vert_temp = 0
        mesh_vertices[i] = numpy.zeros(3)
        for j in range (0, 4):
            mesh_vertices[i] += ( nodes_pos[tet[j]] ) * mesh_itrpltn[i][j]
            vert_temp += T[tet[j]] * mesh_itrpltn[i][j]
        verts_color[i] = to_tmp_color(vert_temp)
    N_f = len(mesh_faces)

def getMeshDisplayInfo():
    global gl_vertices, gl_normals, verts_num, gl_colors, N_f, mesh_vertices
    display_norms = numpy.zeros([N_f * 3, 3])
    display_verts = numpy.zeros([N_f * 3, 3])
    display_colrs = numpy.zeros([N_f * 3, 4])
    for i in range(0, N_f):
        f = mesh_faces[i]
        v1 = mesh_vertices[f[2]] - mesh_vertices[f[1]]
        v2 = mesh_vertices[f[0]] - mesh_vertices[f[1]]
        normal = numpy.cross(v1, v2)
        normal = normal / numpy.linalg.norm(normal)
        for j in range (0, 3):
            display_norms[i * 3 + j] = normal
            display_verts[i * 3 + j] = mesh_vertices[f[j]]
            display_colrs[i * 3 + j] = verts_color[f[j]]
    display_verts_list = display_verts.flatten().tolist()
    display_norms_list = display_norms.flatten().tolist()
    display_colrs_list = display_colrs.flatten().tolist()
    verts_num = len(display_verts_list)
    gl_vertices = (GLfloat * len(display_verts_list))(*display_verts_list)
#     print("display_norms_list: ", display_norms_list)
    gl_normals = (GLfloat * len(display_norms_list))(*display_norms_list)
    gl_colors = (GLfloat * len(display_colrs_list))(*display_colrs_list) 
   
# red being hottest, blue being coldest, green being in the middle
def to_tmp_color(tmp):
    span = T_env - T_init
    middle = T_init + span / 2.0
    red = numpy.array([1., 0, 0])
    green = numpy.array([0, 1., 0])
    blue = numpy.array([0, 0, 1.])
    color = numpy.zeros(3)
    if tmp >= T_env:
        color = red
    elif tmp >= middle:
        ratio = ((tmp - middle) / (span / 2.0)) 
        color = (ratio * red + (1.0 - ratio) * green)
    elif tmp >= T_init:
        ratio = ((tmp - T_init) / (span / 2.0)) 
        color = (ratio * green + (1.0 - ratio) * blue)
    else:
        color = blue
    return numpy.append(color, 1.0)

# different color if node is boundary node
def to_bndry_color(node_index):
    if node_index in bndry_nodes:
        return numpy.array([0.3, 0, 0, 0.5])
    return numpy.array([1., 1., 1., 0.5])

@window.event
def on_draw():
    global frame_num, verts_num

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(4, 5, 5, 0, 0, 0, 0, 1, 0)
    glRotatef(ry, 0, 1, 0)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60., 1., .1, 1000.)
    
    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT)
    
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)

    glEnableClientState(GL_COLOR_ARRAY)

#     the grid
    prepareVerts()
    glVertexPointer(3, GL_FLOAT, 0, gl_vertices)
    glNormalPointer(GL_FLOAT, 0, gl_normals)
    glColorPointer(4, GL_FLOAT, 0, gl_colors)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
    glDrawArrays(GL_TRIANGLES, 0, int(N_f) * 3)
    
    # the mesh
    getMeshVertices()
    getMeshDisplayInfo()
    glVertexPointer(3, GL_FLOAT, 0, gl_vertices)
    glNormalPointer(GL_FLOAT, 0, gl_normals)
    glColorPointer(4, GL_FLOAT, 0, gl_colors)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    glDrawArrays(GL_TRIANGLES, 0, int(N_f) * 3)
    
    glPopClientAttrib()
    if save_imgs:
        mgr = pyglet.image.get_buffer_manager()
        mgr.get_color_buffer().save('images/screenshot'+str(frame_num)+'.png')
        frame_num += 1
    
setup()
# while(True):
#     animate(0)
#     window.dispatch_event('on_draw')
pyglet.clock.schedule(animate)
pyglet.app.run()


# In[ ]:





# In[ ]:




