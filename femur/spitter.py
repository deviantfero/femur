def init_file(filename):
    with open(filename, "w+") as file:
        file.write("Gid Post Results File 1.1\n\n")
        file.write("# encoding utf-8\n\n")
        file.write("GaussPoints \"tet4_element_gp\" ElemType Tetrahedra\n"
                   "Number Of Gauss Points: 4\n"
                   "Natural Coordinates: Given\n"
                   "0.58541    0.138197    0.138197\n"
                   "0.138197    0.58541    0.138197\n"
                   "0.138197    0.138197    0.58541\n"
                   "0.138197    0.138197    0.138197\n"
                   "End gausspoints\n\n")


def get_header(time, what):
    if what is "velocity":
        return (
            "Result \"PRESSURE\" \"Femur\" {} Scalar OnNodes\n"
            "ComponentNames \"X-VELOCITY\", \"Y-VELOCITY\","
            " \"Z-VELOCITY\", \"|VELOCITY|\"\n"
            "Values\n".format(time)
        )
    if what is "pressure":
        return (
            "Result \"PRESSURE\" \"Femur\" {} Scalar OnNodes\n"
            "ComponentNames \"PRESSURE\"\n"
            "Values\n".format(time)
        )


def write_delta(filename, data):
    with open(filename, "a") as file:
        file.write(get_header(data[0], "velocity"))
        for row in data[1]:
            file.write(" {} {} {} {} {}\n".format(*row))
        file.write("End values\n\n")

        file.write(get_header(data[0], "pressure"))
        for row in data[2]:
            file.write(" {} {}\n".format(*row))
        file.write("End values\n\n")
