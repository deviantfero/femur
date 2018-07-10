import sympy as sym
import numpy as npy
from . import matrices_nav_stokes as matrices
from .symbols import *
from sympy.interactive import printing

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

# # force section
# B = velocity * (matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
#                                   matrices.gradNxX)  # * matrices.velocity_arr
# C = (-matrices.Jdet / rho) * (matrices.NT_integrated_respecto_to_epsilon *
#                               matrices.gradNxX)  # * matrices.overlapped_pressure_arr
# D = matrices.gradNxX_gradNxXt #* matrices.velocity_arr
# E = (-matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
#                         matrices.gradNxX)  # * matrices.velocity_arr

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

    #6 times cause i have 6 variables
    assembly_left_mat = sym.Matrix(conn_size * 6, conn_size * 6, ([0] * ((conn_size * 6)**2)))
    variables_mat = variables_vector(conn_size)
    assembly_right_mat = sym.Matrix(conn_size * 6, 1, [0] * conn_size * 6)

    #output nodes
    output_nodes = [x for x in data["nodes"] if x.is_output]
    
    #input nodes
    input_nodes = [x for x in data["nodes"] if x.is_input]     

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

        #assembling left side
        for i in range(len(conn)):
            for j in range(len(conn)):
                assembly_left_mat[conn[i] * 6 : (conn[i] * 6) + 6, 
                                  conn[j] * 6 : (conn[j] * 6) + 6] += local_left_mat[i * 6 : i * 6 + 6, 
                                                                                     j * 6 : j * 6 + 6] 
                assembly_right_mat[conn[i] * 6 : (conn[i] * 6) + 6, : ] += local_right_mat[i * 6 : (i * 6) + 6, : ] 

    #neumann condition array
    neumann_position_array = []
    for i in range(len(data["nodes"])):
        if data["nodes"][i].is_output:
            neumann_position_array += [i] 
   
    #dirchlet condition array
    dirchlet_position_array = []
    for i in range(len(data["nodes"])):
        if data["nodes"][i].is_input:
            dirchlet_position_array += [i]

    #applying dirchlet condition
    velocity_0 = 15

    #adding the velocity to right mat
    for i in range(len(dirchlet_position_array)):
        #adding to right matrix the vx and vy for each node in the entry face of the model
        assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_position_array[i] * 6] * (-velocity_0)
        assembly_right_mat[ : , 0] += assembly_left_mat[ : , dirchlet_position_array[i] * 6 + 1] * (-velocity_0)

    #deleting ecuations solved
    dirchlet_position_array = dirchlet_position_array[::-1]    

    rem = sym.Symbol("rem")
    #removing vx and vy from nodes in entry point
    for i in range(len(dirchlet_position_array)):
        assembly_left_mat[dirchlet_position_array[i] * 6, : ] = sym.Matrix(1, assembly_left_mat.shape[1], [rem] * assembly_left_mat.shape[1])
        assembly_left_mat[dirchlet_position_array[i] * 6 + 1, : ] = sym.Matrix(1, assembly_left_mat.shape[1], [rem] * assembly_left_mat.shape[1])

        assembly_left_mat[ : , dirchlet_position_array[i] * 6] = sym.Matrix(assembly_left_mat.shape[0], 1, [rem] * assembly_left_mat.shape[0])
        assembly_left_mat[ : , dirchlet_position_array[i] * 6] = sym.Matrix(assembly_left_mat.shape[0], 1, [rem] * assembly_left_mat.shape[0])

        variables_mat[dirchlet_position_array[i] * 6, : ] = sym.Matrix(1, 1, [rem])
        variables_mat[dirchlet_position_array[i] * 6 + 1, : ] = sym.Matrix(1, 1, [rem])

        assembly_right_mat[dirchlet_position_array[i] * 6, : ] = sym.Matrix(1, 1, [rem])
        assembly_right_mat[dirchlet_position_array[i] * 6 + 1, : ] = sym.Matrix(1, 1, [rem])

    print('K Matrix shape: \n', assembly_left_mat.shape)
    print('Vars Matrix shape: \n', variables_mat.shape)
    print('B Matrix shape: \n', assembly_right_mat.shape)
