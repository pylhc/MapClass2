# Optimisations to pytpsa

## Profiling

Creating a Map of the FFS, using 'fastcomposition', m1 * m2,
these are the methods that take longer:

      ncalls  tottime  percall  cumtime  percall filename:lineno(function)
       22623   29.481    0.001   39.118    0.002 pytpsa/pol.py:251(fmulpol)
    22518992    6.526    0.000    6.526    0.000 {zip}
       23739    4.261    0.000    6.664    0.000 pytpsa/pol.py:221(addpol)
    31872861    3.501    0.000    3.501    0.000 {method 'get' of 'dict' objects}
      109715    1.762    0.000    3.090    0.000 pytpsa/pol.py:161(truncate)
       27127    0.924    0.000   40.435    0.001 pytpsa/pol.py:235(mulpol)
     2108962    0.589    0.000    0.699    0.000 pytpsa/pol.py:10(abs)
     4684159    0.551    0.000    0.551    0.000 {sum}
       35571    0.206    0.000    0.917    0.000 pytpsa/pol.py:205(mulcoef)
         405    0.194    0.000   47.095    0.116 pytpsa/polmap.py:250(compose)

Where the important column is **tottime** for the total time spent on that method.

## Riccardo de Maria suggestions

### First email

One potential line that can be speed up is:

    tuple([l+m for l,m in zip(i,j)]) in fmulpol

The only line of attack I see is to redo the lowlevel stuff in C, in
particular the arithmetic operations. For instance the line above takes
1us for 6 long exponents, while adding 6 integers should take 20ns, so
there is a factor 50 to gain there.

To do that I estimate two months of work, otherwise it is pointless.In detail you need to write the following functions:

    int exp2i(char[] exp, char [] norder, char [] nvar);
    i2exp(int i, char [] norder, char [] nvar);

which translates the position in the array of coefficient into exponents.

In python, this functionality is implemented in dictionary of tuples,
but in C you could use a generator for i2exp, which fix the ordering,
and a memoized b+tree or hash tables for exp2i.

### Second email

Take this code:

    def ncomp(n, k):
       """Number of composition of n in k
       (n + k - 1)! / n! (k-1)!
       (n+1)  (n + k -1) / 1 ... k
       """
       result = 1
       for i in range(1, k):
           result = result * (n+i) / i
       return result

    def nallcomp(n,k):
       """Sum of the number of composition of i in k for i<=n"""
       return ncomp(n,k+1)

    def comp(n,k):
      if k>1:
        for j in range(n+1):
           for rest in comp(j,k-1):
             yield [n-j]+rest
      else:
        yield [n]

    def allcomp(n,k):
      for nn in range(n+1):
        for exp in comp(nn,k):
           yield exp

    order,nvariables=5,6
    i2exp=list(allcomp(order,nvariables))
    len(i2exp)==nallcomp(order,nvariables)

This code generates all exponents of order<=5 and 6 variables, ordered by
order. You can invert i2exp using a b+tree with nvariable long leaves or an
hash table. Once you have this you can put the coefficient in an array and the
location of the product of two monomial i, j  is just

    exp2i(i2exp(i)+i2exp(j))

### Other suggestions

Start to see how many times and how long it takes to 1) add
polynomial, 2) multiply polynomial. Try to see how much you gain between
m1(m2) and m1*m2 if m1 and m2 are polmap.

I had in mind a modification to add a property called 'binding' to a
polynomial.

Now if you have polynomials with a different signature ("vars") then
they are merged together in new polynomials with a larger "vars". This
is what allows to build the polynomials with normal operations. In
certain cases it would be nice that the coefficient of the polynomials
are polynomials themselves.

The idea is then to specify a property called binding. For any binary
operations like \_\_add\_\_, \_\_mul\_\_ ... if the other argument is a not a
polynomial or has a different binding, then it is treated as a scalar
and you call addcoeff or mulcoef, otherwise you cal addpol and mulpol.

You can also add another method faddpol that you call when you are sure
that  other.vars==self.vars like with mulpol and fmulpol. It should
speed up things a bit.
