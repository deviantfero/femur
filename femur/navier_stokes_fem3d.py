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
# A = matrices.Jdet * matrices.NT_integrated_respecto_to_epsilon * matrices.force_arr
# B = velocity * (matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
#                                   matrices.gradNxX)  # * matrices.velocity_arr
# C = (-matrices.Jdet / rho) * (matrices.NT_integrated_respecto_to_epsilon *
#                               matrices.gradNxX)  # * matrices.overlapped_pressure_arr
# D = matrices.gradNxX_gradNxXt #* matrices.velocity_arr
# E = (-matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
#                         matrices.gradNxX)  # * matrices.velocity_arr

def variables_vector(cant):
    strf = ""
    for i in range(cant):
        strv += "vx{} vy{} vz{} ".format(x)
        strp += "p{} p{} p{} ".format(x)
        strf += "{} {} ".format(strv, strp)

    return sym.symbols(strf)

def navier_stokes(data, force_arr):
    conn_size = len(data["nodes"])

    #6 times cause i have 6 variables
    assembly_mat = sym.Matrix(conn_size * 6, conn_size * 6, ([0] * ((conn_size * 6)**2)))

    #solving each local mat
    for conn in data["connections"]:
        nodes_arr = [data["nodes"][conn[0]].get_position(), 
                     data["nodes"][conn[1]].get_position(), 
                     data["nodes"][conn[2]].get_position(), 
                     data["nodes"][conn[3]].get_position()]

        local_mat = navier_stokes_local(matrices.local_lambda, 
                                              [1, 2, 3], 
                                              data["velocity"], 
                                              data["density"], 
                                              nodes_arr)

        #removing non inexact values (nans) from local solutions (in case determinant of gradXxE is zero)
        where_are_NaNs = npy.isnan(local_mat)
        local_mat[where_are_NaNs] = 0

        for i in range(len(conn)):
            for j in range(len(conn)):
                assembly_mat[conn[i] : conn[i] + 6, conn[j] : conn[j] + 6] += local_mat[i : i + 6, j : j + 6] 

        # print(assembly_mat[0,0])

def navier_stokes_local(
        eq,
        force_arr,
        adv_velocity,
        density,
        nodes_arr):

    # nodes
    xs = [x[0] for x in nodes_arr]
    ys = [x[1] for x in nodes_arr]
    zs = [x[2] for x in nodes_arr]

    result = eq(*xs, *ys, *zs,
              *force_arr,
              min(xs), max(xs),
              min(ys), max(ys),
              min(zs), max(zs),
              density, adv_velocity)

    return result
