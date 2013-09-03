import math

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

def F(a,b):
    """
    New!!! Evaluate primitiveF in zeros and add deltas
    a = sqrt(K)*L
    b = sqrt(K)* (l*)
    """
    def primitiveF(a,b):
        """
        This function calculates the F primitive parameter in Oide effect
        a is the integration limit {sqrt(K)*L, if it has no zeros}
        b is the parameter multiplying cos(phi) {sqrt(K)* (l*)}
        """
        
        a2 = a*a
        a3 = a2*a
        a4 = a3*a
        a5 = a4*a
        a6 = a5*a
        a7 = a6*a

        b2 = b*b
        b3 = b2*b
        b4 = b3*b
        b5 = b4*b
        b6 = b5*b
        b7 = b6*b
        
        sina = math.sin(a)
        sin3a = math.sin(3*a)
        sin5a = math.sin(5*a)
        sin7a = math.sin(7*a)
        cosa = math.cos(a)
        cos3a = math.cos(3*a)
        cos5a = math.cos(5*a)
        cos7a = math.cos(7*a)
      
        sinb = math.sin(b)
        sin3b = math.sin(3*b)
        sin5b = math.sin(5*b)
        sin7b = math.sin(7*b)
        cosb = math.cos(b)
        cos3b = math.cos(3*b)
        cos5b = math.cos(5*b)
        cos7b = math.cos(7*b)

        F = (1/1209600.0
             * (1323*cos5a-675*cos7a
                +a*(378000*sina+21000*sin3a-7560*sin5a)
                +b*(23625*sina+4725*sin3a-14175*sin5a+4725*sin7a
                    +a *(-37800*cos5a)
                    +a2*(-75600*sin3a+226800*sina)
                    )
                +b2*(-49707*cos5a+14175*cos7a
                      +a *(1587600*sina-172200*sin3a+68040*sin5a)
                      )
                +b3*(-80325*sina-144725*sin3a+82215*sin5a-23625*sin7a
                      +a *(37800*cos5a)
                      +a2*(680400*sina-126000*sin3a)
                      )
                +b4*(68985*cos5a-23625*cos7a
                     +a *(2041200*sina-205800*sin3a+37800*sin5a)
                     )
                +b5*(-458325*sina-43225*sin3a-25893*sin5a+14175*sin7a
                      +a *(68040*cos5a)
                      +a2*(680400*sina-25200*sin3a)
                      )
                +b6*(-945*cos5a+4725*cos7a
                      +a *(831600*sina-12600*sin3a-37800*sin5a)
                      )
                +b7*((-354375)*sina + 5425*sin3a-1323*sin5a-675*sin7a
                     +a* (-7560*cos5a)
                     +a2*(226800*sina + 25200*sin3a)
                     )
                +cosa*(b2+1)*4725*(
                    b5*a*80
                    +b4*(155-48*a2)
                    +b3*a*64
                    +b2*(182-96*a2)
                    +b *(-16*a)
                    +a2*(-48)
                    +75
                    )
                +cos3a*(-175)*(
                    +b7*(a*(120))
                     +b6*3*(144*a2+71)
                     +b5*a*744
                     +b4*(720*a2+347)
                     +b3*(a*(-24))
                     +b2*(144*a2-473)
                     +b*(-648)*a
                     +a2*(-144)-31
                     )
                )
             )
        return F

    roots = []
    roots.append(primitiveF(0,b))
    phi = math.atan(-b)+math.pi # starts with n*pi, n=1
    while phi < a:
      roots.append(primitiveF(phi,b))
      phi += math.pi
    roots.append(primitiveF(a,b))
    total = 0
    for i in range(1,len(roots)):
      total += abs(roots[i]-roots[i-1])
    return total
