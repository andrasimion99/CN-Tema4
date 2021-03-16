import math


def read_tridiagonal_matrix_file(read_file):
    try:
        f = open(read_file, "r")
        n = int(f.readline())
        p = int(f.readline())
        q = int(f.readline())

        f.readline()
        a = [0 for _ in range(n)]
        for i in range(n):
            a[i] = float(f.readline())

        f.readline()
        c = [0 for _ in range(n - p)]
        for i in range(len(c)):
            c[i] = float(f.readline())

        f.readline()
        b = [0 for _ in range(n - q)]
        for i in range(len(b)):
            b[i] = float(f.readline())

        return a, b, c
    except OSError:
        print("Error while reading data!")


def read_f(read_file):
    try:
        f = open(read_file, "r")
        n = int(f.readline())

        f.readline()
        a = [0 for _ in range(n)]
        for i in range(n):
            a[i] = float(f.readline())

        return a
    except OSError:
        print("Error while reading data!")


def check_matrix_diag(a, epsilon):
    ok = True

    for elem in a:
        if abs(elem) <= epsilon:
            ok = False
    return ok


def gauss_seidel_method(a, b, c, f, epsilon):
    k_max = 10000
    delta_max = pow(10, 8)
    n = len(a)
    p = n - len(c)
    q = n - len(b)
    # xp = [0 for _ in range(n)]
    # xc = [0 for _ in range(n)]
    xgs = [0 for _ in range(n)]
    k = 0
    delta_x = 0
    while True:
        print("iteratia:", k)
        # xp = xc
        xgs, delta_x = calculate_xgs(a, b, c, f, n, q, p, xgs)
        print(delta_x)
        # delta_x = calculate_deltax(xc, xp)
        k += 1
        if k > k_max or delta_x < epsilon or delta_x > delta_max:
            break
    if delta_x < epsilon:
        print("Am gasit solutie:", xgs)
        return xgs
    else:
        exit("divergenta")


def calculate_xgs(a, b, c, f, n, q, p, xgs):
    norm = 0
    for i in range(n):
        s1 = 0
        for j in range(0, i):
            if j == q + i and i < j:
                # am element in b[i]
                s1 += xgs[j] * b[i]
            elif i - p == j and i > j:
                # am element in c[i-p]
                s1 += xgs[j] * c[i - p]
        s2 = 0
        for j in range(i + 1, n):
            if j == q + i and i < j:
                # am element in b[i]
                s2 += xgs[j] * b[i]
            elif i - p == j and i > j:
                # am element in c[i-p]
                s2 += xgs[j] * c[i - p]
        new_val = (f[i] - s1 - s2) / a[i]
        # print(new_val, xgs[i])
        norm += (new_val - xgs[i]) ** 2
        xgs[i] = new_val
    return xgs, math.sqrt(norm)


def calculate_deltax(xc, xp):
    s = 0
    for i in range(len(xc)):
        s += (xc[i] - xp[i]) ** 2
    return math.sqrt(s)


def calculate_norm(a, b, c, xgs, f):
    n = len(a)
    p = n - len(c)
    q = n - len(b)
    # s = [0 for _ in range(n)]
    norm = 0
    for i in range(n):
        sum = 0
        for j in range(n):
            if j == q + i and i < j:
                # am element in b[i]
                # print(b[i])
                sum += b[i] * xgs[j]
            elif i - p == j and i > j:
                # am element in c[i-p]
                # print(c[i - p])
                sum += c[i - p] * xgs[j]
            elif i == j:
                # print(a[i])
                sum += a[i] * xgs[j]
        norm += (sum - f[i]) ** 2
    return math.sqrt(norm)


if __name__ == '__main__':
    a, b, c = read_tridiagonal_matrix_file("a4.txt")
    f = read_f("f4.txt")
    epsilon = pow(10, -10)
    if not check_matrix_diag(a, epsilon):
        exit("Diagonala principala are elemente nule")

    xgs = gauss_seidel_method(a, b, c, f, epsilon)
    norm = calculate_norm(a, b, c, xgs, f)
    print("norma ||A*xgs -f|| =", norm)
