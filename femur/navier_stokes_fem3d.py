import sympy as sym
import math
import numpy as npy
from . import matrices_nav_stokes as matrices
from .symbols import *
from sympy.interactive import printing

# # force section
# B = velocity * (matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
#                                   matrices.gradNxX)  # * matrices.velocity_arr
# C = (-matrices.Jdet / rho) * (matrices.NT_integrated_respecto_to_epsilon *
#                               matrices.gradNxX)  # * matrices.overlapped_pressure_arr
# D = matrices.gradNxX_gradNxXt #* matrices.velocity_arr
# E = (-matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
#                         matrices.gradNxX)  # * matrices.velocity_arr

def is_on_plane(point, plane_eq):
    return 0 - round(plane_eq.subs([(x, point.x), (y, point.y), (z, point.z)])) == round(0.0001, 3)

def surface_from_three_points(point1, point2, point3):
    vect_p1_p = [x - point1.x, y - point1.y, z - point1.z]
    vect_p1_p2 = [point2.x - point1.x, point2.y - point1.y, point2.z - point1.z]
    vect_p1_p3 = [point3.x - point1.x, point3.y - point1.y, point3.z - point1.z]

    surface_eq = sym.Matrix([vect_p1_p, vect_p1_p2, vect_p1_p3])

    return surface_eq.det()


def gradient(num, dom):
    result = sym.Matrix([])

    for row in dom:
        aux = []
        for col in num:
            aux = aux + [sym.diff(col, row)]
        result = result.row_insert(result.shape[0], sym.Matrix([aux]))

    return result

def overlap_matrix_col(matrix, times):
    result = sym.Matrix([[]])

    for col in range(matrix.shape[1]):
        for row in range(times):
            result = result.col_insert(result.shape[1], matrix.col(col))

    return result

def overlap_matrix_row(matrix, times):
    result = sym.Matrix([[]])

    for row in range(matrix.shape[0]):
        for aux in range(times):
            result = result.row_insert(result.shape[0], matrix.row(row))

    return result

def variables_vector(cant):
    aux = []

    for x in range(cant):
        strv = "vx{} vy{} vz{} ".format(x, x, x)
        strp = "p{} p{} p{} ".format(x, x, x)
        strf = "{} {} ".format(strv, strp)
        aux += list(sym.symbols(strf))

    return sym.Matrix(aux)

def navier_stokes_local(
        left,
        right,
        force_arr,
        adv_velocity,
        density,
        nodes_arr):

    # nodes
    xs = [x[0] for x in nodes_arr]
    ys = [x[1] for x in nodes_arr]
    zs = [x[2] for x in nodes_arr]

    result_left = left(*xs, *ys, *zs,
                       min(xs), max(xs),
                       min(ys), max(ys),
                       min(zs), max(zs),
                       density, adv_velocity)

    result_right = right(*xs, *ys, *zs, *force_arr)


    return result_left, result_right

def navier_stokes(data, force_arr):
    conn_size = len(data["nodes"])

    #6 times because i have 6 variables
    assembly_left_mat = sym.Matrix(conn_size * 6, conn_size * 6, ([0] * ((conn_size * 6)**2)))
    variables_mat = variables_vector(conn_size)
    assembly_right_mat = sym.Matrix(conn_size * 6, 1, [0] * conn_size * 6)

    #output nodes
    output_nodes = [x.id for x in data["nodes"] if x.is_output]
    
    #input nodes
    input_nodes = [x.id for x in data["nodes"] if x.is_input]     

    #solving each local mat
    for conn in data["connections"]:
        nodes_arr = [data["nodes"][conn[0]].get_position(), 
                     data["nodes"][conn[1]].get_position(), 
                     data["nodes"][conn[2]].get_position(), 
                     data["nodes"][conn[3]].get_position()]
 
        local_left_mat, local_right_mat = navier_stokes_local(matrices.left_side_lambda,
                                                              matrices.right_side_lambda,
                                                              [1, 2, 3], 
                                                              data["velocity"], 
                                                              data["density"], 
                                                              nodes_arr)

        # #removing non inexact values (nans) from local solutions (in case determinant of gradXxE is zero)
        local_left_mat = npy.nan_to_num(local_left_mat)
        local_right_mat = npy.nan_to_num(local_right_mat)

        #assembling left and right side
        for i in range(len(conn)):
            for j in range(len(conn)):
                assembly_left_mat[conn[i] * 6 : (conn[i] * 6) + 6, 
                                  conn[j] * 6 : (conn[j] * 6) + 6] += local_left_mat[i * 6 : i * 6 + 6, 
                                                                                     j * 6 : j * 6 + 6] 
                assembly_right_mat[conn[i] * 6 : (conn[i] * 6) + 6, : ] += local_right_mat[i * 6 : (i * 6) + 6, : ] 

    #applying no slip conditions
    conns_mat = data["connections"]
    noslip_elems = data["noslip_elems"]
    noslip_connections = [data["connections"][x] for x in data["noslip_elems"]]
    noslip_nodes = [x for x in data["nodes"] if x.is_noslip]

    #finding plane and removing nodes that belong to plane (only for no slip nodes)
    count = 0
    dirchlet_extra_nodes = [] #velocity = 0 with this ones
    while len(noslip_elems) > 0:
        curr_elem = noslip_elems[0]
        nodes_pos_in_elem = conns_mat[curr_elem]

        #getting the node objects in an array
        nodes_vals_in_elem = []
        for i in nodes_pos_in_elem:
            nodes_vals_in_elem += [data["nodes"][i]]

        #noslip nodes for this iteration (with this i create the plane)
        curr_noslip_nodes = [x for x in nodes_vals_in_elem if x.is_noslip]
        plane = sym.Poly(surface_from_three_points(*curr_noslip_nodes[:3]))

        nodes_in_plane = [x.id for x in data["nodes"] if is_on_plane(x, plane)]

        # #removing all nodes in plane
        # normal_vector = sym.Matrix([plane.coeffs()[:3]])
        # if normal_vector.shape[1] < 2:
        #     normal_vector.col_insert(1, sym.Matrix([0]))
        #     normal_vector.col_insert(2, sym.Matrix([0]))
        # else: 
        #     if normal_vector.shape[1] < 3:
        #         normal_vector.col_insert(2, sym.Matrix([0]))
        #
        # unit_vector = normal_vector / math.sqrt(normal_vector[0, 0]**2 + normal_vector[0, 1]**2 + normal_vector[0, 2]**2)
        # noslip_velocity_vector = unit_vector * 638

        # #removing all no_slip elements that contain this nodes
        # elem_match = []
        #
        # for i, elem in enumerate(noslip_connections):
        #     for node in elem:
        #         cond = 0
        #         for nip in nodes_in_plane:
        #             cond += data["nodes"][node] == nip
        #     elem_match.append(cond > 0)
        #
        # count += 1
        #
        # new_elems = [noslip_elems[x] for x in range(len(noslip_elems)) if not elem_match[x]]
        dirchlet_extra_nodes += nodes_in_plane
        noslip_elems = noslip_elems[1:]


    dirchlet_extra_nodes = list(set(dirchlet_extra_nodes))

    #dirchlet condition array
    dirchlet_position_array = []
    for i in range(len(data["nodes"])):
        if data["nodes"][i].is_input:
            dirchlet_position_array += [i]
    dirchlet_position_array.sort()

    #applying dirchlet condition
    velocity_0 = -638

    #adding the initial velocity to right mat
    for i in range(len(dirchlet_position_array)):
        #adding to right matrix the vx and vy for each node in the entry face of the model
        assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_position_array[i] * 6] * (-velocity_0)
        assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_position_array[i] * 6 + 1] * (-velocity_0)

    # for i in range(len(dirchlet_extra_nodes)):
    #     assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_extra_nodes[i] * 6] * (-velocity_0) #vx
    #     assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_extra_nodes[i] * 6 + 1] * (-velocity_0) #vy
    #     assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_extra_nodes[i] * 6 + 1 + 1] * (-velocity_0) #vz

    
    #deleting ecuations solved
    dirchlet_extra_nodes.sort()
    dirchlet_position_array += dirchlet_extra_nodes
    dirchlet_position_array += input_nodes
    dirchlet_position_array = list(set(dirchlet_position_array))
    dirchlet_position_array.sort()
    dirchlet_position_array = dirchlet_position_array[::-1]    

    rem = sym.Symbol("rem")
    #removing vx and vy from nodes in entry point
    for i in range(len(dirchlet_position_array)):
        assembly_left_mat[dirchlet_position_array[i] * 6, : ] = sym.Matrix(1, assembly_left_mat.shape[1], [rem] * assembly_left_mat.shape[1])
        assembly_left_mat[dirchlet_position_array[i] * 6 + 1, : ] = sym.Matrix(1, assembly_left_mat.shape[1], [rem] * assembly_left_mat.shape[1])

        assembly_left_mat[ : , dirchlet_position_array[i] * 6] = sym.Matrix(assembly_left_mat.shape[0], 1, [rem] * assembly_left_mat.shape[0])
        assembly_left_mat[ : , dirchlet_position_array[i] * 6 + 1] = sym.Matrix(assembly_left_mat.shape[0], 1, [rem] * assembly_left_mat.shape[0])

        variables_mat[dirchlet_position_array[i] * 6, : ] = sym.Matrix(1, 1, [rem])
        variables_mat[dirchlet_position_array[i] * 6 + 1, : ] = sym.Matrix(1, 1, [rem])

        assembly_right_mat[dirchlet_position_array[i] * 6, : ] = sym.Matrix(1, 1, [rem])
        assembly_right_mat[dirchlet_position_array[i] * 6 + 1, : ] = sym.Matrix(1, 1, [rem])


    #removing rows and cols from K matrix
    m, n = assembly_left_mat.shape
    rows = [i for i in range(m) if any(assembly_left_mat[i, j] != rem for j in range(n))]
    cols = [j for j in range(n) if any(assembly_left_mat[i, j] != rem for i in range(m))]

    assembly_left_mat = assembly_left_mat[rows, cols]

    #removing rows from variables
    m, n = variables_mat.shape
    rows = [i for i in range(m) if any(variables_mat[i, j] != rem for j in range(n))]
    cols = [j for j in range(n) if any(variables_mat[i, j] != rem for i in range(m))]

    variables_mat = variables_mat[rows, cols]

    #removing rows from B matrix
    m, n = assembly_right_mat.shape
    rows = [i for i in range(m) if any(assembly_right_mat[i, j] != rem for j in range(n))]
    cols = [j for j in range(n) if any(assembly_right_mat[i, j] != rem for i in range(m))]
        
    assembly_right_mat = assembly_right_mat[rows, cols]
 
    print('No slip Nodes: \n', len(noslip_nodes))
    print('K Matrix shape: \n', assembly_left_mat.shape)
    print('Vars Matrix shape: \n', variables_mat.shape)
    print('B Matrix shape: \n', assembly_right_mat.shape)

    #solving time 
    delta_time = data["time_delta"]
    end_time = data["end_time"]
    f = sym.lambdify(*variables_mat, variables_mat)
    curr_x = f(*[10 for x in range(variables_mat.shape[0])])
    for i in range(0, end_time, delta_time):
        next_x = (sym.eye(assembly_left_mat.shape[0]) + delta_time * assembly_left_mat).inv() * (curr_x + delta_time * assembly_right_mat)
        print(next_x)
        curr_x = next_x

