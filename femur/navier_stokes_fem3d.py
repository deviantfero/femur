import sympy as sym
from symbols import * 
import matrices_nav_stokes as matrix

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
            result = result.row_insert(result.shape[1], matrix.row(aux))

    return result

# calculate section of navier-stokes equation with force envolved
def force_section(nt, var_arr):
    var1, var2, var3 = var_arr

    # integrating
    fst = sym.integrate(sym.expand(nt), (var1, 0, 1))
    snd = sym.integrate(fst, (var2, 0, 1))
    trd = sym.integrate(snd, (var3, 0, 1))

    return trd

# calculate section of navier-stokes equation with velocity envolved


def velocity_section(nt, gradNxX, var_arr):
    var1, var2, var3 = var_arr

    result = nt * gradNxX

    # integrating
    fst = sym.integrate(result, (var1, 0, 1))
    snd = sym.integrate(fst, (var2, 0, 1))
    trd = sym.integrate(snd, (var3, 0, 1))

    return trd


def pressure_section(nt, gradNxX, var_arr):
    var1, var2, var3 = var_arr

    result = nt * gradNxX

    # integrating
    fst = sym.integrate(sym.expand(result), (var1, 0, 1))
    snd = sym.integrate(fst, (var2, 0, 1))
    trd = sym.integrate(snd, (var3, 0, 1))

    return trd

# calculate section of navier-stokes equation with velocity envolved and
# double gradient


def weak_velocity_section(nt, gradNxX, x_arr):
    var1, var2, var3 = x_arr

    result = gradNxX.T * gradNxX

    # integrating
    fst = sym.integrate(sym.expand(result), (var1, a, b))
    snd = sym.integrate(fst, (var2, c, d))
    trd = sym.integrate(snd, (var3, f, g))
    return trd


def navier_stokes_local(
        force_arr,
        pressure_arr,
        adv_velocity,
        density,
        nodes_arr):

    # nodes
    node1, node2, node3, node4 = nodes_arr

    # x1, y1, z1 = node1
    # x2, y2, z2 = node2
    # x3, y3, z3 = node3
    # x4, y4, z4 = node4

    # arrays
    isoparametric_arr = [e, n, l]
    velocity_arr = sym.Matrix(
        [vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx4, vy4, vz4])

    # approximation of interpolation
    x_interpolated = (x2 - x1) * e + (x3 - x1) * n + (x4 - x1) * l + x1
    y_interpolated = (y2 - y1) * e + (y3 - y1) * n + (y4 - y1) * l + y1
    z_interpolated = (z2 - z1) * e + (z3 - z1) * n + (z4 - z1) * l + z1

    X_arr = [x_interpolated, y_interpolated, z_interpolated]

    # shape functions for isoparametric domain
    N1 = 1 - e - n - l
    N2 = e
    N3 = n
    N4 = l

    N = [N1, N2, N3, N4]

    # arrays for creating the N matrix
    N_1 = [N1, 0, 0]
    N_2 = [N2, 0, 0]
    N_3 = [N3, 0, 0]
    N_4 = [N4, 0, 0]

    # underhand N transymosed matrix
    NT_mat = sym.Matrix([
        N_1, N_1[::-2] + [0], N_1[::-1],
        N_2, N_2[::-2] + [0], N_2[::-1],
        N_3, N_3[::-2] + [0], N_3[::-1],
        N_4, N_4[::-2] + [0], N_4[::-1]])

    # underhand N matrix
    N_mat = NT_mat.T

    # gradient of N respect to epsilon
    gradNxE = gradient(N, isoparametric_arr)
    overlapped_gradNxE = overlap_matrix_col(gradNxE, 3)

    gradXxE = gradient(X_arr, isoparametric_arr)
    overlapped_gradXxE = overlap_matrix_col(gradXxE, 3)

    # gradient of N respect to x
    gradXxE = sym.simplify(gradient(X_arr, isoparametric_arr))
    gradNxX = sym.simplify(gradXxE.inv()) * overlapped_gradNxE

    #NT_mat_integrated_respect_to_epsilon = sym.integrate(sym.integrate(sym.integrate(NT_mat, (e, 0, 1)), (n, 0, 1)), (l, 0, 1))

    return (gradNxX.T * gradNxX) * (b - a) * (c - d) * (g - f)

    # pressure arr
    overlapped_pressure = overlap_matrix_row(sym.Matrix(pressure_arr) / 3, 3)

    # jacobian
    J = gradient(X_arr, isoparametric_arr).T

    # determinant
    Jdet = J.det()


    # force section
    A = Jdet * force_section(NT_mat, isoparametric_arr) * sym.Matrix(force_arr)
    B = (Jdet * adv_velocity) * velocity_section(NT_mat,
                                                 gradNxX, isoparametric_arr) * velocity_arr
    C = ((-Jdet / density) * pressure_section(NT_mat,
                                              gradNxX, isoparametric_arr)) * overlapped_pressure
    D = weak_velocity_section(NT_mat, gradNxX, [x, y, z]) * velocity_arr
    E = (-Jdet) * velocity_section(NT_mat,
                                   gradNxX, isoparametric_arr) * velocity_arr

    return A + B + C + D + E


# testing values
if __name__ == "__main__":
    print(navier_stokes_local([1, 2, 3], [4, 5, 6, 7], 5, 1000, [
          [31, 2, 2], [41, 56, 3], [12, 45, 59], [10, 10, 4]]))
