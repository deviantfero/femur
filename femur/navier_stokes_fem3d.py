import sympy as sym
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

    return eq(*xs, *ys, *zs,
              *force_arr,
              min(xs), max(xs),
              min(ys), max(ys),
              min(zs), max(zs),
              density, adv_velocity)
