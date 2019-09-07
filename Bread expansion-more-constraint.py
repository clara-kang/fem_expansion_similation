#!/usr/bin/env python
# coding: utf-8

# # Expansion with Collision

# In[1]:


import math
import numpy
import pyglet
from pyglet.gl import *

# bread radius
radius = 4 # in grid unit

# grid cell size
grid_size = 10

# grid size
X = 10
Y = 10

# center of bread circle
center = numpy.array([X/2, Y/2])

# triangles by node index
triangles = numpy.zeros([0, 3], dtype=int)

# grid 1 for node exist
grid = numpy.zeros([X, Y], dtype=int)

# node
nodes = numpy.zeros([0, 2])

# map grid to node index
grid2i = numpy.zeros([X, Y], dtype=int)

# the index of nodes on the boundary
bndry_nodes = numpy.zeros(0, dtype=int)

# get circle intersection with y and x grid lines
itrsctn_y = numpy.zeros([Y, 2])
itrsctn_x = numpy.zeros([X, 2])
for y in range (0, Y):
    if (y - center[1])**2 > radius ** 2:
        continue
    itrsctn_x1 = -math.sqrt(radius ** 2 - (y - center[1])**2) + center[0]
    itrsctn_x2 = math.sqrt(radius ** 2 - (y - center[1])**2) + center[0]
    itrsctn_y[y] = numpy.array([itrsctn_x1, itrsctn_x2])
for x in range (0, X):
    if (x - center[0])**2 > radius ** 2:
        continue
    itrsctn_y1 = -math.sqrt(radius ** 2 - (x - center[0])**2) + center[1]
    itrsctn_y2 = math.sqrt(radius ** 2 - (x - center[0])**2) + center[1]
    itrsctn_x[x] = numpy.array([itrsctn_y1, itrsctn_y2])
    
# get nodes
cnt = 0
for y in range (0, Y):
    if itrsctn_y[y][0] != 0:
        left_is_whole = False
        right_is_whole = False
        # check the left intersection is fraction
        if itrsctn_y[y][0] % 1 != 0:
            grid[y][math.floor(itrsctn_y[y][0])] = 1
            grid2i[y][math.floor(itrsctn_y[y][0])] = cnt
            bndry_nodes = numpy.append(bndry_nodes, cnt)
            nodes = numpy.append(nodes, [[itrsctn_y[y][0], y]], axis=0)
            cnt += 1
        else:
            left_is_whole = True
        # check the right intersection is fraction
        if itrsctn_y[y][1] % 1 != 0:
            grid[y][math.ceil(itrsctn_y[y][1])]= 1
            grid2i[y][math.ceil(itrsctn_y[y][1])] = cnt
            bndry_nodes = numpy.append(bndry_nodes, cnt)
            nodes = numpy.append(nodes, [[itrsctn_y[y][1], y]], axis=0)
            cnt += 1
        else:
            right_is_whole = True
        # fill the nodes between left and right intersection
        if left_is_whole:
            bndry_nodes = numpy.append(bndry_nodes, cnt)
        for x in range (math.ceil(itrsctn_y[y][0]), math.floor(itrsctn_y[y][1]) + 1):
            grid[y][x] = 1
            grid2i[y][x] = cnt
            nodes = numpy.append(nodes, [[x, y]], axis=0)
            if x == math.floor(itrsctn_y[y][1]) and right_is_whole and math.ceil(itrsctn_y[y][0]) != math.floor(itrsctn_y[y][1]):
                bndry_nodes = numpy.append(bndry_nodes, cnt)
            cnt += 1

# get nodes on the border for each x line
for x in range (0, X):
    if itrsctn_x[x][0] % 1 != 0:
        grid[math.floor(itrsctn_x[x][0])][x] = 1
        grid2i[math.floor(itrsctn_x[x][0])][x] = cnt
        bndry_nodes = numpy.append(bndry_nodes, cnt)
        nodes = numpy.append(nodes, [[x, itrsctn_x[x][0]]], axis=0)
        cnt += 1
    if itrsctn_x[x][0] % 1 != 0:
        grid[math.ceil(itrsctn_x[x][1])][x] = 1
        grid2i[math.ceil(itrsctn_x[x][1])][x] = cnt
        bndry_nodes = numpy.append(bndry_nodes, cnt)
        nodes = numpy.append(nodes, [[x, itrsctn_x[x][1]]], axis=0)
        cnt += 1
        
# construct triangle elements       
for x in range (0, X-1):
    for y in range (1, Y):
        if grid[y][x] != 0:
            if grid[y][x+1] != 0 and grid[y-1][x] != 0: # upper right triangle
                triangles = numpy.append(triangles, [[grid2i[y][x], grid2i[y][x+1], grid2i[y-1][x]]], axis=0)
            elif grid[y][x+1] != 0 and grid[y-1][x+1] != 0:
                triangles = numpy.append(triangles, [[grid2i[y][x], grid2i[y][x+1], grid2i[y-1][x+1]]], axis=0)
for x in range (1, X):
    for y in range (0, Y-1):
        if grid[y][x] != 0:
            if grid[y][x-1] != 0 and grid[y+1][x] != 0: # upper left triangle
                triangles = numpy.append(triangles, [[grid2i[y][x], grid2i[y][x-1], grid2i[y+1][x]]], axis=0)
            elif grid[y][x-1] != 0 and grid[y+1][x-1] != 0:
                triangles = numpy.append(triangles, [[grid2i[y][x], grid2i[y][x-1], grid2i[y+1][x-1]]], axis=0)

print("bndry_nodes: ", bndry_nodes)

numpy.set_printoptions(precision=1)

# number of nodes
N_n = nodes.shape[0]

# number of triangles
N_t = triangles.shape[0]

# damping constant
damp_constant = 3

# the displacements
u = numpy.zeros([N_n, 2])

# the velocities
v = numpy.zeros([N_n, 2])

# the force on each node
F = numpy.zeros([N_n, 2])

# B matrices
B = numpy.zeros([N_t, 3, 6])

# B transpose mult by D
B_T_D = numpy.zeros([N_t, 6, 3])

# stiffness matrix
K = numpy.zeros([int(N_n*2), int(N_n*2)])
 
# gravity
g = numpy.array([0, -0.98])

# mass of each node
m = 1 #kg

# external force
F_ex = numpy.zeros([N_n, 2])

# yield force magnitude
y = 3

# initialize D-matrix for Hooke's Law
D = numpy.identity(3)
D[0][1] = 0.5
D[1][0] = 0.5
D *= 7.4e1

# calculate B matrix for triangles
def updateBMatrices():
    global B        
    M = numpy.zeros([3,3]) # for calculating determinant d
    b, c, x, y = numpy.zeros(3), numpy.zeros(3), numpy.zeros(3), numpy.zeros(3) # shape function constants
    for i in range (0, N_t):
        triangle = triangles[i]
        for j in range (0, 3):
            x[j] = nodes[triangle[j]][0] # easier to access triangle's vertices
            y[j] = nodes[triangle[j]][1]
        for j in range (0, 3):
            M[j] = numpy.array([1, x[j], y[j]])
        d = numpy.linalg.det(M)
        b[0] = y[1] - y[2]
        b[1] = y[2] - y[0]
        b[2] = y[0] - y[1]
        c[0] = x[2] - x[1]
        c[1] = x[0] - x[2] 
        c[2] = x[1] - x[0]
        b /= d
        c /= d
        B[i][0] = numpy.array([b[0], 0, b[1], 0, b[2], 0])
        B[i][1] = numpy.array([0, c[0], 0, c[1], 0, c[2]])
        B[i][2] = 0.5 * numpy.array([c[0], b[0], c[1], b[1], c[2], b[2]])
  
# assemble stiffness matrix
def updateK():
    global K, B
    for i in range (0, N_t):
        B_T_D[i] = numpy.matmul(B[i].T, D)
        K_l = numpy.matmul(B_T_D[i], B[i])
        for j in range (0, 3):
            for k in range (0, 3):
                K[2*triangles[i][j]][2*triangles[i][k]] += K_l[2*j][2*k]
                K[2*triangles[i][j]+1][2*triangles[i][k]] += K_l[2*j+1][2*k]
                K[2*triangles[i][j]][2*triangles[i][k]+1] += K_l[2*j][2*k+1]
                K[2*triangles[i][j]+1][2*triangles[i][k]+1] += K_l[2*j+1][2*k+1]
    
updateBMatrices()
updateK()

to_kelvin = 273.15

# conductivity of material
k = 2.36e-2

# strain-temperature constant
alpha = 1e-2

# heat capacitance of material
c = 0.9

# heat radiation constant
ems_sig = 2.5e1

# initial temperature
T_init = 0 # celcius

# env temperature
T_env = 150 # celcius

# each node's temperature
T = T_init * numpy.ones(N_n)
print("T: ", T)

# force due to heat
h = numpy.zeros([N_n, 2])

# each node's heat flux
q_ex = numpy.zeros(N_n)

# stiffness for heat for each element
K_e = numpy.zeros([N_t, 3, 3])

# calculate B matrix for triangles
def init_Ke():
    global K_e
    M = numpy.zeros([3,3]) # for calculating determinant d
    N = numpy.zeros([3,3]) # shape functions
    b, c, x, y = numpy.zeros(3), numpy.zeros(3), numpy.zeros(3), numpy.zeros(3) # shape function constants
    for i in range (0, N_t):
        triangle = triangles[i]
        for j in range (0, 3):
            x[j] = nodes[triangle[j]][0] # easier to access triangle's vertices
            y[j] = nodes[triangle[j]][1]
        for j in range (0, 3):
            M[j] = numpy.array([1, x[j], y[j]])
        A = numpy.linalg.det(M) / 2.0
        b[0] = y[1] - y[2]
        b[1] = y[2] - y[0]
        b[2] = y[0] - y[1]
        c[0] = x[2] - x[1]
        c[1] = x[0] - x[2] 
        c[2] = x[1] - x[0]
        K_l = numpy.array([b, c])
        K_e[i] = (k / (4.0 * A ** 2.0)) * numpy.matmul(K_l.T, K_l)

init_Ke()

# assemble stiffness matrix
K_h = numpy.zeros([N_n, N_n])
for i in range (0, N_t):
    for j in range (0, 3):
        for k in range (0, 3):
            K_h[triangles[i][j]][triangles[i][k]] += K_e[i][j][k]
print ("K_h: \n", K)

def get_bndry_force():
    global T
    b_force = numpy.zeros(N_n)
    for i in bndry_nodes:
        b_force[i] = ems_sig * ( T_env - T[i] )
    return b_force


# calculate "heat" force
def updateh():
    global K, B, h, B_T_D
    h.fill(0)
    for i in range (0, N_t):
        avg_T = ( T[triangles[i][0]] + T[triangles[i][1]] + T[triangles[i][2]] ) / 3.0
        strain_h = numpy.array([alpha * avg_T, alpha * avg_T, 0])
        h_l = numpy.matmul(B_T_D[i], strain_h)
        for j in range (0, 3):
            h[triangles[i][j]][0] += h_l[j*2]
            h[triangles[i][j]][1] += h_l[j*2+1]
    
# floor y position
floor_y = 0.5

# left wall position
lwall_x = - 0.5 

# right wall position
rwall_x = 9 

# height of wall
wall_h = 10

# time step, don't use actual dt which is not stable
dt = 0.01

# apply impact
def col_dtctn():
    global F_ex, v, u
    F_ex.fill(0) # reset external force
    F_ex += g
    nodes_pos = nodes + u
    for i in range(0, N_n):
        node = nodes[i]
        pos = nodes_pos[i]
        if pos[1] <= floor_y:
            u[i][1] = floor_y - node[1]
            v[i][1] = - 0.5 * v[i][1]
        if pos[0] <= lwall_x and pos[1] < wall_h:
            u[i][0] = lwall_x - node[0]
            v[i][0] = - 0.5 * v[i][0]
        elif pos[0] >= rwall_x  and pos[1] < wall_h:
            u[i][0] = rwall_x - node[0]
            v[i][0] = - 0.5 * v[i][0]
        
# red being hottest, blue being coldest, green being in the middle
def to_tmp_color(tmp):
    span = T_env - T_init
    middle = T_init + span / 2.0
    red = numpy.array([255, 0, 0])
    green = numpy.array([0, 255, 0])
    blue = numpy.array([0, 0, 255])
    if tmp >= T_env:
        return red.astype(int)
    elif tmp >= middle:
        ratio = ((tmp - middle) / (span / 2.0)) 
        return (ratio * red + (1.0 - ratio) * green).astype(int)
    elif tmp >= T_init:
        ratio = ((tmp - T_init) / (span / 2.0)) 
        return (ratio * green + (1.0 - ratio) * blue).astype(int)
    else:
        return blue.astype(int)


def animate(dt_actual):
    global u, v, F, K, dt, F_ex, nodes
    global K_h, T, c
    
    F.fill(0) # reset force
    updateh();
    def getF(v, u):
        global F_ex, K, N_n
        F_in = numpy.array( (numpy.matmul(K, u.flatten())).reshape([N_n,2])) # calculate internal force
        F = numpy.zeros([N_n, 2])
        F += F_ex # add external force
        F -= v * damp_constant # account for damping
        F += h # add heat force
        F -= F_in # minus internal force
        print ("F: ", F)
        return F
#     # start of Heun's method
#     F_e = getF(v, u)
#     v_e = v + dt * F_e / m
#     u_e = u + dt * v    

#     v += dt * (F_e + getF(v_e, u_e)) / (2.0*m)
#     u += dt * (v + v_e) / 2.0
#     # end of Heun's method
    v += dt * getF(v, u) / m # semi-implicit euler
    u += dt * v # semi-implicit euler
    col_dtctn() # do collision detection to get degrees of freedom
    
    q_in = -numpy.matmul(K_h, T) # internal heat transfer
    q_ex = get_bndry_force() # external heat transfer
    q = q_in + q_ex
    dQ = q * dt
    T += (1.0 / (c * m)) * dQ 
    
#     print("T: ", T)

# create window
window = pyglet.window.Window(500, 500)
pyglet.gl.glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
pyglet.gl.glLineWidth(2)

def drawTriangles():
    # make it look bigger
    origin = numpy.array([200, 30])
    nodes_to_draw = origin + grid_size * (nodes + u)
    # get colro based on temperature
    nodes_color = (numpy.zeros([N_n, 3])).astype(int)
    for i in range(0, N_n):
        nodes_color[i] = to_tmp_color(T[i])
    # draw
    pyglet.graphics.draw_indexed(nodes_to_draw.shape[0], pyglet.gl.GL_TRIANGLES,
    triangles.astype(int).flatten(),
    ('v2f', tuple(nodes_to_draw.flatten()) ),
    ('c3B', tuple(nodes_color.flatten()))
)
    
@window.event
def on_draw():
    window.clear()
    drawTriangles()
    
u.fill(0)
v.fill(0)
T.fill(0) 
pyglet.clock.schedule_interval(animate, .01)
pyglet.app.run()   


# In[ ]:





# In[ ]:




