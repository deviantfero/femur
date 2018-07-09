import fileinput
from .node import Node


def parse_info():
    """
    Parses GiD generated information, places
    return: dictionary with GiD information
    the dictionary has the following keys
    :connections: a connection table for all elements
    :nodes: a list of Node objects
    :node_count: amount of nodes in model
    :time_delta: the step of time in simulation
    :velocity: initial velocity
    """
    problem_info = {"connections": []}
    nodes = []
    stage = 0

    for line in fileinput.input():
        if ";" in line:
            stage += 1
        elif "#" not in line:
            line = line.replace("\n", "").split(" ")
            data = [float(x) for x in line if x is not ""]
            if data and stage == 0:
                problem_info["node_count"] = data[0]
                problem_info["time_delta"] = data[2]
                problem_info["velocity"] = data[3]
                problem_info["density"] = data[4]
            if data and stage == 1:
                nodes.append(Node(*data[1:]))
                if len(nodes) == problem_info["node_count"]:
                    problem_info["nodes"] = nodes
            if data and stage == 2:
                problem_info["connections"] += [[int(x - 1) for x in data[1:]]]
            if data and stage == 3:
                problem_info["nodes"][int(data[0] - 1)].is_noslip = True
            if data and stage == 4:
                problem_info["nodes"][int(data[0] - 1)].is_input = True
            if data and stage == 5:
                problem_info["nodes"][int(data[0] - 1)].is_output = True
    return problem_info
