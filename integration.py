def simpson(f, a, b, n):
    """
    Simpsons rule

    :param function f: a unary function
    :param int a: start point
    :param int b: end point
    :param int n: no. integration intervals (must be even)

    :return: result of the integral
    """

    h = (float(b) - a) / n
    S = f(float(a))

    for i in range(1, n, 2):
        x = a + h * i
        S += 4 * f(x)

    for i in range(2, n-1, 2):
        x = a + h * i
        S += 2 * f(x)

    S += f(float(b))
    F = h * S / 3

    return F


def trap(f, a, b, n):
    """
    Trapezium rule

    :param function f: a unary function
    :param int a: start point
    :param int b: end point
    :param int n: no. integration intervals (must be even)

    :return: result of the integral
    """

    h = (float(b) - float(a)) / float(n)
    S = f(float(a)) + f(float(b))

    for i in range(1, n):
        x = a + h * i
        S += 2 * f(x)

    F = h * S / 2

    return F
