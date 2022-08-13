from math import sqrt


def isarealnumber(a):
    try:
        a = float(a)
        if str(a) == "nan" or str(a) == "inf" or str(a) == "-inf":
            return False
        return True
    except (ValueError):
        return False


def eq_str_editor(s):
    if not isinstance(s, str):
        raise ValueError
    s = s.replace("1(", "(")
    s = s.replace("^1", "")
    s = s.replace("(x^0)", "")
    s = s.replace("(x)", "x", 1)
    s = s.rstrip(" ")
    s = s.replace(" ", " + ")
    s = s.replace(" + -", " - ")
    return s


class Point:
    def __init__(self, x, y):
        if not (isarealnumber(x) and isarealnumber(y)): raise ValueError
        self.x = x
        self.y = y
        self.coordinates = (self.x, self.y)
        self.PO = self.dToB((0, 0))

    def __str__(self):
        return "(x, y) = ({0}, {1})".format(self.x, self.y)

    def __repr__(self):
        return self.__str__()

    def dToB(self, B):
        if isinstance(B, Point):
            return sqrt((self.x - B.x) ** 2 + (self.y - B.y) ** 2)
        elif isinstance(B, tuple):
            x_B = B[0]
            y_B = B[1]
            return sqrt((self.x - x_B) ** 2 + (self.y - y_B) ** 2)
        else:
            return float("NaN")


class Polynom():
    #  here we assume a polynom a*x^n + b*x^(n-1) + c*x^(n-2) ..... = 0 is about to be being processed.
    #  so a, b, c ..... coefficients should be given by the correct order, where all coefficients are real numbers.
    #  Important: Zeroes must be included. Such as, for a polynom x^3, the given arguments must be Polynom(1, 0, 0, 0).
    def __init__(self, *coefficients):
        try:
            for ce in coefficients:
                if isinstance(ce, float) or isinstance(ce, int):
                    continue
                else:
                    raise ValueError
            self.coefficients = coefficients
            self.d = len(self.coefficients) - 1
            for coef in self.coefficients:
                if coef == 0:
                    self.d -= 1
                else:
                    break
            self.degree = self.d
            self.y0 = self.coefficients[-1]

        except ValueError:
            print("Coefficients must be real numbers.")

    def __str__(self):
        return self.info().replace("The polynom is", "A polynom with the equation", 1)

    def __repr__(self):
        return self.__str__()

    def f(self, x):
        if x == 0:
            return self.coefficients[-1]
        if x == 1:
            return sum(self.coefficients)
        sum = 0
        temp_degree = len(self.coefficients) - 1
        for c in self.coefficients:
            sum += c * (x ** temp_degree)
            temp_degree -= 1
        return sum

    def info(self):
        info = ""
        t_d = len(self.coefficients) - 1
        for c in self.coefficients:
            if c == 0:
                t_d -= 1
            else:
                if t_d != 0:
                    info += "{0}(x^{1}) ".format(c, t_d)
                    t_d -= 1
                else:
                    if c != 0:
                        info += str(c)
        return "This polynom is P(x) = " + eq_str_editor(info)


class Parabola(Polynom):
    def __init__(self, a1, a2, a3):
        if a1 == 0:
            raise ZeroDivisionError
        super().__init__(a1, a2, a3)
        self.a1 = a1
        self.delta = a2 ** 2 - (4 * a1 * a3)
        self.r = (-a2 / (2 * a1))
        self.k = self.f(self.r)
        self.vertex = Point(self.r, self.k)

    def __str__(self):
        return self.info().replace("polynom", "parabola", 1)

    def __repr__(self):
        return self.__str__()

    def the_roots(self):
        if self.delta < 0:
            return (float('nan'), float('nan'))
        elif self.delta == 0:
            return (self.k, self.k)
        else:
            return ((self.k + (sqrt(self.delta) / (2 * self.a1))), (self.k - (sqrt(self.delta) / (2 * self.a1))))


class Cubic(Polynom):
    def __init__(self, p, q, r, t):
        if p == 0:
            raise ZeroDivisionError
        super().__init__(p, q, r, t)
        self.derivative_of_cubic = Parabola(3 * p, 2 * q, r)
        if self.derivative_of_cubic.delta < 0:
            self.vertices = (None, None)
        elif self.derivative_of_cubic.delta == 0:
            temp_x_of_vertex = self.derivative_of_cubic.the_roots()[0]
            self.v_1 = Point(temp_x_of_vertex, self.f(temp_x_of_vertex))
            self.vertices = (self.v_1, self.v_1)
        else:
            temp_x1_of_v1 = self.derivative_of_cubic.the_roots()[0]
            temp_x2_of_v2 = self.derivative_of_cubic.the_roots()[1]
            v1 = Point(temp_x1_of_v1, self.f(temp_x1_of_v1))
            v2 = Point(temp_x2_of_v2, self.f(temp_x2_of_v2))
            self.vertices = (v1, v2)
        self.delta_0 = q ** 2 - (3 * p * r)
        self.delta_1 = 2 * (q ** 3) - 9 * (p * q * r) + 27 * (p ** 2) * t
        self.p, self.q, self.r, self.t = p, q, r, t
        self.delta = (4 * (self.delta_0 ** 3) - (self.delta_0 ** 2)) / (27 * p * p)

    def __str__(self):
        return self.info().replace("polynom", "cubic", 1)

    def __repr__(self):
        return self.__str__()

    def the_roots(self):
        if self.delta_1 == 0 and self.delta_0 == 0:
            c_r = (-self.q / (3 * self.p))
            return (c_r, c_r, c_r)


class Line(Polynom):
    #  Equation type is the type of equation of the line that will be given to the init function.
    #  L: Linear Equation, such as ax + by + c = 0, a, b anc c are being expected. This type is set by default.
    #  Expected tuple for L mode is (a, b, c), where a,b and c are real numbers.
    #  (if b is given as zero, it will be treated specifically, not as a polynomial.)
    #  P: Polynomial, such as y = mx + n, where m is the slope of the line and n is the interception of line with y axis.
    #  Expected tuple for P mode is (m, n), where m and n are real numbers.
    #  !!!_a, b, c, m, or n, all values should be given in the tuple regardless of whether they are zero or not!_!!!
    def __init__(self, t_of_coefs, equation_type="L"):
        if equation_type == "L":
            self.a, self.b, self.c = t_of_coefs

            if self.a == self.b == 0:
                raise ValueError("Both the x's and y's coefficients cannot be zero simultaneously.")
            if self.a == 0:
                self.m = 0
                self.n = - self.c / self.b
                self.type = "H"  # meaning that this line is 'H'orizontal.
                super().__init__(self.n)
            elif self.b == 0:
                self.m = float("inf")
                self.n = None
                self.degree = None
                self.y0 = None
                self.type = "V"  # Similarly, meaning that this line is 'V'ertical.

            else:
                self.m = - self.a / self.b
                self.n = - self.c / self.b
                self.type = str(self.m)
                super().__init__(self.m, self.n)
        elif equation_type == "P":
            self.m, self.n = t_of_coefs
            self.a, self.b, self.c = self.m, -1, self.n

            super().__init__(self.m, self.n)
        else:
            raise ValueError("Given argument for \'equation_type\':[0} is unexpected.".format(equation_type))
        if self.m == 0:
            self.type = "H"
        elif self.m == float("inf"):
            self.type = "V"
        else:
            self.type = str(self.m)

    def f(self, x):
        if self.m == 0:
            return self.n
        elif self.m == float("inf"):
            return None

        else:
            return (x * self.m) + self.n

    def info(self):
        _info = "This is a horizontal line whose equation is y = "
        if isarealnumber(self.m):
            if self.m == 0:
                return _info + str(self.n)
            elif self.m == 1:
                _info = _info.replace(" horizontal", "") + str("x +" + str(self.n)).replace("+ -", " - ")
                return _info
            elif self.m == -1:
                return str("This is a line whose equation is y = -x +" + str(self.n)).replace("+ -", " - ")
            else:
                return "This is a line whose equation is y = {0}x + {1}".format(self.m, self.n).replace("+ -", "- ")
        else:
            return "This is a vertical line such that x = " + str(-self.c / self.a)

    def __str__(self):
        return self.info()

    def __repr__(self):
        return self.__str__()

    def wrt_l2(self, l2):
        if not isinstance(l2, Line):
            raise ValueError("Given argument must be a line.")
        if self.type == l2.type:
            if self.f(1) == l2.f(1):
                return "C"  # meaning that these lines are coinciding.
            else:
                return "P"  # meaning that these lines are only parallel but not coincident.
        else:
            return "I"  # meaning that these lines are intercepting.

    def interception_with(self, l2):
        if not isinstance(l2, Line):
            raise ValueError("Given argument must be another line.")
        _s = self.wrt_l2(l2)
        if _s == "C" or _s == "P":
            return ((None, None), _s)
        else:
            if self.a != 0:
                y = (((self.c * l2.a) / self.a) - l2.c) / (l2.b - ((l2.a * self.b) / self.a))
                x = -((self.b * y) + self.c) / self.a
                temp_p = Point(x, y)
                return temp_p
            else:
                y = -self.c / self.b
                x = -(l2.c + l2.b * y) / l2.a
                temp_p = Point(x, y)
                return temp_p

    def inv_f(self, x):
        if self.m == float("inf") or self.m == 0:
            return None
        else:
            return (x - self.n) / self.m
