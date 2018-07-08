import sympy as sym
import matrices_nav_stokes as matrices
from symbols import *
from sympy.interactive import printing

# get the gradient of num respect to dom


def gradient(num, dom):
    result = sym.Matrix([])

    for row in dom:
        aux = []
        for col in num:
            aux = aux + [sym.diff(col, row)]
        result = result.row_insert(result.shape[0], sym.Matrix([aux]))

    return result

# overlap matrix's col


def overlap_matrix_col(matrix, times):
    result = sym.Matrix([[]])

    for col in range(matrix.shape[1]):
        for row in range(times):
            result = result.col_insert(result.shape[1], matrix.col(col))

    return result

# overlap matrix's row


def overlap_matrix_row(matrix, times):
    result = sym.Matrix([[]])

    for row in range(matrix.shape[0]):
        for aux in range(times):
            result = result.row_insert(result.shape[0], matrix.row(row))

    return result


# force section
A = matrices.Jdet * matrices.NT_integrated_respecto_to_epsilon  # * matrices.force_arr
B = velocity * (matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
                                  matrices.gradNxX)  # * matrices.velocity_arr
C = (-matrices.Jdet / rho) * (matrices.NT_integrated_respecto_to_epsilon *
                              matrices.gradNxX)  # * matrices.overlapped_pressure_arr
D = matrices.gradNxX_gradNxXt * matrices.velocity_arr
E = (-matrices.Jdet) * (matrices.NT_integrated_respecto_to_epsilon *
                        matrices.gradNxX)  # * matrices.velocity_arr


def navier_stokes_local(
        eq,
        force_arr,
        pressure_arr,
        adv_velocity,
        density,
        nodes_arr):

    # nodes
    xs = [x[0] for x in nodes_arr]
    ys = [x[1] for x in nodes_arr]
    zs = [x[2] for x in nodes_arr]

    return eq(*xs, *ys, *zs,
              *pressure_arr,
              min(xs), max(xs),
              min(ys), max(ys),
              min(zs), max(zs),
              *force_arr, density, velocity)

print(matrices.local_mat.shape)

# testing values
if __name__ == "__main__":
    for i in range(1000):
        print(navier_stokes_local(matrices.local_lambda, [1, 2, 3], [4, 5, 6, 7], 5, 1000, [
            [31, 2, 2], [41, 56, 3], [12, 45, 59], [10, 10, 4]]).shape)
