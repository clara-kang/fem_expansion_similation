#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy

N = 4

numpy.set_printoptions(precision=2, threshold=numpy.inf)

# map grid position to node index
grid2node = numpy.zeros([N+1, N+1, N+1], dtype=int)

# nodes position
nodes = numpy.zeros([(N+1)**3, 3])

# nodes in the left half:
left_half_nodes = numpy.zeros(0)

#nodes on the boundary
bndry_nodes = numpy.zeros(0, dtype=int)

# tetrahedra elements, each containing nodes index
tetrahedra = numpy.zeros([N*N*N*5, 4], dtype=int)

# triangles on the boundary/surface, each triangle's nodes index
srfce_trgls = numpy.zeros([N*N*6*2, 3], dtype=int)

def createNodes():
    global bndry_nodes, nodes, grid2node, left_half_nodes
    cnt = 0
    for i in range (0, N+1):
        for j in range (0, N+1):
            for k in range (0, N+1):
                nodes[cnt] = numpy.array([i, j, k])
                grid2node[i][j][k] = cnt
                if i <= (N+1)/2:
                    left_half_nodes = numpy.append(left_half_nodes, cnt)
                if i == 0 or j == 0 or k == 0 or i == N or j == N or k == N:
                    bndry_nodes = numpy.append(bndry_nodes, cnt)
                cnt += 1

def createTets():
    cnt = 0
    for i in range (0, N):
        for j in range (0, N):
            for k in range (0, N):
                origin = grid2node[i][j][k]
                back = grid2node[i][j][k+1]
                right = grid2node[i+1][j][k]
                back_right = grid2node[i+1][j][k+1]
                up = grid2node[i][j+1][k]
                up_back = grid2node[i][j+1][k+1]
                up_right = grid2node[i+1][j+1][k]
                up_right_back = grid2node[i+1][j+1][k+1]
                
                tetrahedra[cnt] = numpy.array([origin, back_right, back, up_back])
                tetrahedra[cnt + 1] = numpy.array([right, back_right, origin, up_right])
                tetrahedra[cnt + 2] = numpy.array([up_back, up_right, up, origin])
                tetrahedra[cnt + 3] = numpy.array([up_right_back, up_right, up_back, back_right])
                tetrahedra[cnt + 4] = numpy.array([up_back, up_right, back_right, origin])
                cnt += 5
                
def createSurfTriangles():
    cnt = 0
    v1 = v2 = v3 = v4 = 0
    def store(k):
        nonlocal cnt, v1, v2, v3, v4
        if k == 0:
            srfce_trgls[cnt] = numpy.array([v1, v2, v3])
            srfce_trgls[cnt+1] = numpy.array([v1, v3, v4])
        elif k == N:
            srfce_trgls[cnt] = numpy.array([v1, v3, v2])
            srfce_trgls[cnt+1] = numpy.array([v1, v4, v3])
        cnt += 2
    for i in range (0, N):
        for j in range (0, N):
            for k in [0, N]:
                # bot and top
                v1 = grid2node[i][k][j]
                v2 = grid2node[i+1][k][j]
                v3 = grid2node[i+1][k][j+1]
                v4 = grid2node[i][k][j+1]
                store(k)
                # left and right
                v1 = grid2node[k][i][j]
                v2 = grid2node[k][i][j+1]
                v3 = grid2node[k][i+1][j+1]
                v4 = grid2node[k][i+1][j]
                store(k)
                # front and back
                v1 = grid2node[i][j][k]
                v2 = grid2node[i][j+1][k]
                v3 = grid2node[i+1][j+1][k]
                v4 = grid2node[i+1][j][k]
                store(k)

createNodes()
createTets()
createSurfTriangles()


# # The physics

# In[2]:


# number of nodes
N_n = nodes.shape[0]

# number of triangles
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
g = numpy.array([0, -0.98, 0])

# mass of each node
m = 1 #kg

# external force
F_ex = numpy.zeros([N_n, 3])

# modulus of elasticity
E = 1e2

# poisson ratio
p_r = 0.3

# floor y position
floor_y = 0

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
    global coeffs
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
            u[i][1] = floor_y - node[1]
            v[i][1] = - 0.5 * v[i][1] # apply some damping

def updateForce(dt_actual):
    global u, v, F, K, dt, F_ex, h
    F.fill(0) # reset force
    F_in = numpy.array( (numpy.matmul(K_s, u.flatten())).reshape([N_n,3])) # calculate internal force
    F += F_ex # add external force
    F -= v * damp_constant # account for damping
    F += h # add heat force
    F -= F_in # minus internal force
    v += dt * F / m
    u += dt * v
#     print ("u: ", u)

getCoeffs()  
updateBsMatrices()
updateKs()

# def moveDown(dt_actual):
#     global u
#     u += numpy.array([0, -0.1, 0])
# moveDown(0)


# In[3]:


# conductivity of material
k_cndct = 2.36e-3

# strain-temperature constant
alpha = 1e-3

# heat capacitance of material
c = 0.9

# heat radiation constant
ems_sig = 2.5e-1

# initial temperature
T_init = 150 #celcius

# env temperature
T_env = 200 # celcius

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
    T_test = 100 * numpy.ones(N_n)
    T_test[0] *= 2.0
    res = numpy.matmul(K_h, T_test)
    print("res: ", res)
                
# calculate heat exchange on boundary nodes              
def getBndryFlux():
    global T
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

# In[4]:


from pyglet.gl import *
import pyglet

gl_vertices = list()
gl_normals = list()
gl_colors = list()

# number of vertices
verts_num = 0

# render surface only
render_surf = False

# render half of cube
render_half = True

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
    
    glLightfv(GL_LIGHT0, GL_POSITION, vec(5, 5, 10))
    glLightfv(GL_LIGHT0, GL_SPECULAR, vec(1, 1, 1, 1))
    glLightfv(GL_LIGHT0, GL_DIFFUSE, vec(1, 1, 1, 1))
    
    glLightfv(GL_LIGHT1, GL_POSITION, vec(0, 0, 5))
    glLightfv(GL_LIGHT1, GL_DIFFUSE, vec(.5, .5, .5, 1))
    glLightfv(GL_LIGHT1, GL_SPECULAR, vec(1, 1, 1, 1))
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, vec(193./255., 217./255., 1.0, 0.5))
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, vec(1, 1, 1, 0.5))
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50)

def prepareVerts():
    global gl_vertices, gl_normals, gl_colors, verts_num
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
                colors.extend(to_tmp_color(T[tetrahedra[j][indices[i][k]]]))
    verts_num = len(vertices)
    gl_vertices = (GLfloat * len(vertices))(*vertices)
    gl_normals = (GLfloat * len(normals))(*normals) 
    gl_colors = (GLfloat * len(colors))(*colors) 
    
def prepareSurfVerts():
    global gl_vertices, gl_normals, verts_num
    normals = []
    vertices = []
    nodes_pos = nodes + u
    for j in range (0, srfce_trgls.shape[0]):
            n = numpy.cross(nodes_pos[srfce_trgls[j][2]] - nodes_pos[srfce_trgls[j][1]], 
                            nodes_pos[srfce_trgls[j][0]] - nodes_pos[srfce_trgls[j][1]])
            n = n / numpy.linalg.norm(n)
            normals.extend(n.tolist() * 3)
            for k in range (0, 3):
                vertices.extend(nodes_pos[srfce_trgls[j][k]])
    verts_num = len(vertices) 
    gl_vertices = (GLfloat * len(vertices))(*vertices)
    gl_normals = (GLfloat * len(normals))(*normals)   


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
    return numpy.append(color, 0.5)

    
# amound to translate to align cube center with origin
def toCenter():
    return [-float(N)/2.] * 3

@window.event
def on_draw():
    global frame_num, verts_num
    if render_surf:
        prepareSurfVerts()
    else:
        prepareVerts()
#     glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) # line mode
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(8, 10, 10, 0, 0, 0, 0, 1, 0)
    glRotatef(ry, 0, 1, 0)
    glTranslatef(*toCenter())
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60., 1., .1, 1000.)
    
    glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT)
    
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)

    glVertexPointer(3, GL_FLOAT, 0, gl_vertices)
    glNormalPointer(GL_FLOAT, 0, gl_normals)

    if not render_surf:
        glEnableClientState(GL_COLOR_ARRAY)
        glColorPointer(4, GL_FLOAT, 0, gl_colors)

    glDrawArrays(GL_TRIANGLES, 0, int(verts_num / 3))
    
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




