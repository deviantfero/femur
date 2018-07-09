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

    for conn in data["connections"]:

        nodes_arr = [data["nodes"][conn[0]].get_position(), 
                     data["nodes"][conn[1]].get_position(), 
                     data["nodes"][conn[2]].get_position(), 
                     data["nodes"][conn[3]].get_position()]

        print(nsfem.navier_stokes_local(matrices.local_lambda, [1, 2, 3], 5, 1000, nodes_arr))


if __name__ == "__main__":
    main()
