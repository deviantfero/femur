import numpy as npy
import sympy as sym
from . import matrices_nav_stokes as matrices
from . import navier_stokes_fem3d as nsfem
from . import parser


def main():
    data = parser.parse_info()

    # print('data-nc', data["node_count"])
    # print('data-td', data["time_delta"])
    # print('data-v', data["velocity"])
    # print('data-d', data["density"])
    # print('connection', data["connections"])

    nsfem.navier_stokes(data, [1, 2, 3])

if __name__ == "__main__":
    main()
