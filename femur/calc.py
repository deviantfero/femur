import numpy as np


def local_fem(k, q, xs, ys):
    a, b = [np.min(xs), np.max(xs)]
    c, d = [np.min(ys), np.max(ys)]

    left = np.array([[ys[1] - ys[2], xs[2] - xs[1]],
                     [ys[2] - ys[0], xs[0] - xs[2]],
                     [ys[0] - ys[1], xs[0] - xs[1]]])
    # print("izquierda:\n", left)
    # print("izquia-trans:\n", np.transpose(left))
    left = left.dot(np.transpose(left))

    ldet = np.linalg.det(np.array([[xs[1] - xs[0], ys[1] - ys[0]],
                                   [xs[2] - xs[0], ys[2] - ys[0]]]))**2

    res = k * (b-a) * (d-c) / ldet
    rdet = np.linalg.det(np.array([[xs[1] - xs[0], xs[2] - xs[0]],
                                   [ys[1] - ys[0], ys[2] - ys[0]]]))
    right = q * rdet * np.array([[0], [0.5], [0.5]])

    # print("ldet:", ldet, "left:\n", left)
    # print("rdet:", rdet)
    return [res * left, right]
