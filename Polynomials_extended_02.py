from math import cos, sin, pi, isnan
from random import uniform

# Links to all mathematical backround used here which are Aberth's method to find all roots simultaneously and Lagrange's improved bounds on all roots of a polynomial
# https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Bounds_on_all_roots
# https://en.wikipedia.org/wiki/Aberth_method
# I tested this code with some polynomials that I created by myself with their roots on paper. It does seem to work!

class Polynom:

    # an*x^n + an-1*x^n-1 + ... a3x^3 + a2x^2 + a1^x + a0 = P(x) is the polynom should be given. P(x) can be defined
    # as *var_name* = Polynom(an, an-1, ..., a3, a2, a1, a0) in the given order, the initial zeros will be eliminated,
    # which can and will decrease the order of the polynomial. E.g. P = Polynom(0, 0, 0, 1, 2, 3) will be regarded as
    # P(x) = x^2 + 2x + 3

    def __init__(self, *coefficients):  # coefficients = [an, an-1, an-2, ..., a3, a2, a1, a0]
        self.a = list(coefficients)
        index = 0
        while index == 0:  # initial zeroes are being eliminated here.
            if self.a[index] == 0:
                del self.a[index]
            else:
                index += 1
        del index
        self.a.reverse()  # creating the main coefficient list in which coefficients are listed as a = [a0, a1, a2, ..., an-1, an]
        self.degree = len(self.a) - 1  # degree of the polynomial is determined from the number of coefficients
        if not bool(
                self.a):  # an exception occurs when the list 'a' is empty, it means P(x) = 0, whose degree is not defined, so interpreted as float 'nan' here.
            self.degree = float('nan')
            self.a = [0]

    #  derivative of the polynomial is being created here
    def derivative_of(self):
        temp_arr = []
        for i in range(1, len(self.a)):
            temp_arr.append(self.a[i] * i)
        temp_arr.reverse()
        temp_pol = Polynom(*temp_arr)
        return temp_pol

    def eval_at(self, x):  # Evaluation of the polynomial at x.
        temp_sum = 0
        for i in range(len(self.a)):
            temp_sum += (self.a[i]) * (x ** i)

        return temp_sum


def cmplx_sqrt(z):  # function to calculate the square-root of floats, integers and complex numbers.
    if isinstance(z, (float, int)):
        if z >= 0:
            return z ** .5
        else:
            return complex(0, abs(z) ** .5)

    else:
        r = abs(z)
        z_ = (r ** 0.5) * ((z + r) / abs(z + r))
        return z_


def get_upper(P):  # Function to get the upper limit of the absolute values of the roots of a polynomial P. I made use
    # of the Lagrange's improved bound calculations for polynomials.
    # Argument: a Polynomial object -> Returns: Upper bound for the polynomial object
    temp_arr, an = [], P.a[-1]
    if P.degree == 0:
        return float('nan')
    for i in range(1, P.degree + 1):
        temp_arr.append(abs(P.a[-(i + 1)] / an) ** (1 / i))
    temp_arr.sort()
    if len(temp_arr) > 1:
        U = temp_arr[-1] + temp_arr[-2]  # Upper Bound
    else:
        U = 2 * temp_arr[0]
    return float(U)


def get_lower(P):  # To find a lower bound for the roots of the polynomial object passed as parameter into the function
    temp_a = P.a.copy()
    reverseP = Polynom(*temp_a)
    rU = get_upper(reverseP)
    if rU == 0:
        return 0
    elif isnan(rU):
        return rU
    else:
        return 1 / rU


#  At this point, one can have an upper and lower bounds for the roots of any polynomial.

def getrandim(n, absU, absL, p=12, onlyReal=False):
    # Using random.uniform and polar coordinates to create pseudo-random complex numbers in an interval.
    # Due to the fact that any z = a + bi = r*cos(theta) + i*r*sin(theta)
    # Here, r is the modulo of z which is created randomly, theta is a randomly generated angle between 0 and 2*pi.
    complex_num_set = []
    if onlyReal:
        return list(x := uniform(absL, absU) for i in range(n))
    while len(complex_num_set) < n:
        rand_r = uniform(absL, absU)
        rand_theta = uniform(0, 2 * pi)
        complex_num_set.append(complex(round(rand_r * cos(rand_theta), p), round(rand_r * sin(rand_theta), p)))

    return complex_num_set


#  A function to calculate the average of the offsets for the approximations of the roots to determine the average error
# It actually just calculates the average mean of all the numbers in the list given as the parameter of the function.
def get_abs_average(_iter):
    temp_sum = 0
    for i in _iter:
        temp_sum += i
    return temp_sum / len(_iter)


#  Using the Aberth's method to calculate the zeros of a given polynomial P. It uses the random set of complex numbers
#  that is created by getrandim and Polynomial's derivative. Its error is 1e-15 by default. The error is calculated
#  from the average mean of the offset numbers.
#  Note: For polynomials of degree 0, 1 and 2; the Aberth's method is not being used. Instead, straightforward algebraic
#  equations are being used. E.g. quadratic formula for polynomials of degree 2.
def get_zeros(P, error=1e-15, prec=3):
    if P.degree == 0:
        return float('nan')
    elif P.degree == 1:
        return tuple(-(P.a[0] / P.a[1]))
    elif P.degree == 2:
        discriminant = P.a[1] ** 2 - 4 * P.a[0] * P.a[2]
        x_1, x_2 = (-P.a[1] + cmplx_sqrt(discriminant)) / (2 * P.a[2]), (-P.a[1] - cmplx_sqrt(discriminant)) / (
                2 * P.a[2])
        return x_1, x_2
    else:
        _z = getrandim(P.degree, get_upper(P), get_lower(P))
        _w = []
        p_ = P.derivative_of()
        err = error + 1
        while abs(err) > error:
            for k in range(P.degree):
                temp_sum = 0
                for j in range(P.degree):
                    if j == k:
                        continue
                    else:
                        temp_sum += 1 / (_z[k] - _z[j])

                ratio = P.eval_at(_z[k]) / p_.eval_at(_z[k])

                _w.append(ratio / (1 - (ratio * temp_sum)))
            err = get_abs_average(_w)

            for index in range(P.degree):
                _z[index] = _z[index] - _w[index]
            _w = []

        for r in range(len(_z)):

            if abs(_z[r].real) < 0.1 * error:
                _z[r] = complex(0, round(_z[r].imag, prec))

            else:
                _z[r] = complex(_z[r].real, round(_z[r].imag, prec))

            if abs(_z[r].imag) < 0.1 * error:
                _z[r] = round(_z[r].real, prec)
            else:
                _z[r] = complex(round(_z[r].real, prec), _z[r].imag)

        _z.sort(key=abs)

        return _z

