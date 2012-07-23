# SIMPSONS RULE (1/3)
# f = function (e.g. I5)
# a = start point, b = end point, n = no. integration intervals

##### n MUST BE EVEN #####

def simpson(f, a, b, n):
    h = (float(b)-a)/n;
    S = f(float(a))

    for i in range(1, n, 2):
        x = a + h*i
        S += 4*f(x)

    for i in range(2, n-1, 2):
        x = a + h*i
        S += 2*f(x)

    S += f(float(b))
    F = h*S/3

    return F

# TRAPEZIUM RULE
# f = function (e.g. I5)
# a = start point, b = end point, n = no. integration intervals
def trap(f, a, b, n):
    h = (float(b)-float(a))/float(n)
    S = f(float(a)) + f(float(b))

    for i in range (1, n):
        x = a + h*i
        S += 2*f(x)

    F = h*S/2

    return F
