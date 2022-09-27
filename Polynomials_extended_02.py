from math import cos, sin, pi, isnan
from random import uniform


# Links to all mathematical background used here which are Aberth's method to find all roots simultaneously and
# Lagrange's improved bounds on all roots of a polynomial
# https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Bounds_on_all_roots
# https://en.wikipedia.org/wiki/Aberth_method
# I tested this code with some polynomials that I created by myself with their roots on paper. It does seem to work!



def cmplx_sqrt(z):  # function to calculate the square-root of floats, integers and complex numbers.
    if isinstance(z, (float, int)):
        if z >= 0:
            return z ** .5
        else:
            return complex(0, (-z) ** .5)

    else:
        r = abs(z)
        z_ = (r ** 0.5) * ((z + r) / abs(z + r))
        return z_


def getrandim(n, absU, absL, p=12, onlyReal=False):
    # Using random.uniform and polar coordinates to create pseudo-random complex numbers in an interval.
    # Due to the fact that any z = a + bi = r*cos(theta) + i*r*sin(theta)
    # Here, r is the modulo of z which is created randomly, theta is a randomly generated angle between 0 and 2*pi.
    if isnan(absU) or isnan(absL):
        return None
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


# Polynomial class to represent and utilise the mathematical objects called polynomials is defined here.
class Polynomial:

    # an*x^n + an-1*x^n-1 + ... a_3x^3 + a_2x^2 + a_1^x + a_0 = P(x) is the polynomial should be given. P(x) can be defined
    # as *var_name* = Polynomial(a_n, a_n-1, ..., a_2, a_1, a_0) in the given order, the initial zeros will be eliminated,
    # which can and will decrease the order of the polynomial. E.g. P = Polynom(0, 0, 0, 1, 2, 3) will be regarded as
    # P(x) = x^2 + 2x + 3

    def __init__(self, *coefficients):  # coefficients = [an, an-1, an-2, ..., a3, a2, a1, a0]

        self.a = list(coefficients)

        while True:
            try:
                if self.a[0] == 0:
                    del self.a[0]
                else:
                    break
            except IndexError:
                break

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
        temp_pol = Polynomial(*temp_arr)
        return temp_pol

    def eval_at(self, x):  # Evaluation of the polynomial at x.
        temp_sum = 0
        for i in range(len(self.a)):
            temp_sum += (self.a[i]) * (x ** i)

        return temp_sum

    def get_upper(self):  # Function to get the upper limit of the absolute values of the roots of a polynomial P. I made use
        # of the Lagrange's improved bound calculations for polynomials.
        # Argument: a Polynomial object -> Returns: Upper bound for the polynomial object
        temp_arr, an = [], self.a[-1]
        if self.degree == 0:
            return float('nan')
        elif isnan(self.degree):
            return float('nan')
        for i in range(1, int(self.degree) + 1):
            temp_arr.append(abs(self.a[-(i + 1)] / an) ** (1 / i))
        temp_arr.sort()
        if len(temp_arr) > 1:
            U = temp_arr[-1] + temp_arr[-2]  # Upper Bound
        else:
            U = 2 * temp_arr[0]
        return float(U)

    def get_lower(self):  # To find a lower bound for the roots of the polynomial object passed as parameter into the function
        temp_a = self.a.copy()
        reverseP = Polynomial(*temp_a)
        rU = reverseP.get_upper()
        if rU == 0:
            return 0
        elif isnan(rU):
            return 0
        else:
            return 1 / rU
        #  At this point, one can have an upper and lower bounds for the roots of any polynomial.

    #  Using the Aberth's method to calculate the zeros of a given polynomial P. It uses the random set of complex numbers
    #  that is created by getrandim and Polynomial's derivative. Its error is 1e-15 by default. The error is calculated
    #  from the average mean of the offset numbers.
    #  Note: For polynomials of degree 0, 1 and 2; the Aberth's method is not being used. Instead, straightforward algebraic
    #  equations are being used. E.g. quadratic formula for polynomials of degree 2.
    def get_zeros(self, error=1e-15, prec=3):
        if self.degree == 0:  # for constant non-zero polynomial, all real line or none is zero. So it returns 'nan.'
            return float('nan')
        elif isnan(self.degree):  # for zero polynomial, all real line is zero, so again, it returns 'nan.'
            return float('nan')
        elif self.degree == 1:  # By the equation 0 = a1x + a0 => x = - a0/a1, where a1 =/= 0.
            return -(self.a[0] / self.a[1])
        elif self.degree == 2:  # basic quadratic formula with the ability of returning complex numbers.
            discriminant = self.a[1] ** 2 - 4 * self.a[0] * self.a[2]
            x_1, x_2 = (-self.a[1] + cmplx_sqrt(discriminant)) / (2 * self.a[2]), (-self.a[1] - cmplx_sqrt(discriminant)) / (
                    2 * self.a[2])
            return x_1, x_2
        elif sum(self.a) == self.a[-1]:  # for the polynomials whose only leading term's coefficient is non-zero.
            return [0]
        else:
            _z = getrandim(self.degree, self.get_upper(), self.get_lower())
            _w = []
            p_ = self.derivative_of()
            err = error + 1
            while abs(err) > error:
                for k in range(int(self.degree)):
                    temp_sum = 0
                    for j in range(int(self.degree)):
                        if j == k:
                            continue
                        else:
                            temp_sum += 1 / (_z[k] - _z[j])

                    ratio = self.eval_at(_z[k]) / p_.eval_at(_z[k])

                    _w.append(ratio / (1 - (ratio * temp_sum)))
                err = get_abs_average(_w)

                for index in range(int(self.degree)):
                    _z[index] = _z[index] - _w[index]
                _w = []
            for r in range(len(_z)):
                # this piece of code gets rid of the approximation(s) (real or imaginary parts of these numbers)
                # that are smaller than the tenth of th error given to the function, default error = 1x10^-15
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

    # compositing the self polynomial and another polynomial Q, in mathematical terms, this function returns R(x) = P(Q(x))
    def comp(self, Q):
        if isnan(self.degree) or self.degree == 0:
            return self
        if isnan(Q.degree):
            return Polynomial(self.a[0])
        elif Q.degree == 0:
            return Polynomial(self.eval_at(Q.a[0]))
        else:
            i = 0
            terms = []
            for a in self.a:
                terms.append(Polynomial(a) * (Q ** i))
                i += 1
            del i
            R = Polynomial(0)
            for _p in terms:
                R += _p
            return R

    #  Common polynomial arithmetics
    def __neg__(self):
        temp_arr = list(term := -self.a[-index] for index in range(1, len(self.a) + 1))
        return Polynomial(*temp_arr)

    def __add__(self, Q):
        temp_a = []
        index = 0
        if isnan(self.degree):
            return Q
        elif isnan(Q.degree):
            return self
        for index in range(int(min(self.degree, Q.degree)) + 1):
            temp_a.append(self.a[index] + Q.a[index])

        if self.degree == Q.degree:
            temp_a.reverse()
            return Polynomial(*temp_a)
        elif self.degree > Q.degree:
            temp_a += self.a[index + 1:]
            temp_a.reverse()
            return Polynomial(*temp_a)
        else:
            temp_a += Q.a[index + 1:]
            temp_a.reverse()
            return Polynomial(*temp_a)

    def __sub__(self, Q):
        return self + (-Q)

    def __mul__(self, Q):
        # scalar multiplication as follows
        if isinstance(Q, (float, int)):
            terms = []
            for a in self.a:
                terms.append(a*Q)
            terms.reverse()
            return Polynomial(*terms)
        if isnan(self.degree) or isnan(Q.degree):
            return Polynomial(0)
        k = self.degree + Q.degree
        temp_q_a = Q.a + [0] * (k - Q.degree)
        temp_p_a = self.a + [0] * (k - self.degree)

        _c = []
        for index in range(k + 1):
            temp_sum = 0
            for j in range(index + 1):
                temp_sum += temp_p_a[j] * temp_q_a[index - j]
            _c.append(temp_sum)
        _c.reverse()
        return Polynomial(*_c)

    def __pow__(self, power):
        if not isinstance(power, int):
            raise TypeError("{0} is not supported as a power of a polynomial.".format(type(power)))
        elif power < 0:
            raise ValueError("Negative integers are not supported as polynomial exponents.")
        if power == 0:
            return Polynomial(1)
        elif power == 1:
            return self
        else:
            # if power >= 2 and is even, then it goes into a recursion and using the formula P^power = (P ^ power/2) * (P ^ power/2)
            if power % 2 == 0:
                return (self ** int(power / 2)) * (self ** int(power / 2))
            # if power is odd, it becomes P * (P ^ (power-1)/2 ) * (P ^ (power-1)/2) )
            else:
                return self * (self ** int((power - 1) / 2)) * (self ** int((power - 1) / 2))

