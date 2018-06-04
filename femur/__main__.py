from . import calc


def main():
    k, q = [float(n) for n in input("input k and q: ").split()]
    xs = [float(x) for x in input("input x values: ").split()]
    ys = [float(y) for y in input("input y values: ").split()]
    print(calc.local_fem(k, q, xs, ys))


if __name__ == "__main__":
    main()
