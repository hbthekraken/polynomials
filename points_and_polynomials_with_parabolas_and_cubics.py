from math import sqrt


def isarealnumber(a):
    try:
        a = float(a)
        if str(a) == "nan" or str(a) == "inf" or str(a) == "-inf":
            return False
        return True
    except (ValueError):
        return False


def str_editor(s):
    if not isinstance(s, str):
        raise ValueError
    s = s.replace("1(", "(")
    s = s.replace("^1", "")
    s = s.replace("(x^0)", "")
    s = s.replace("(x)", "x", 1)
    s = s.rstrip(" ")
    s = s.replace(" ", " + ")
    s = s.replace(" + -", " -")
    return s

def approximation(C, a):
    if not (isinstance(C, Cubic) and isinstance(a, float)):
        raise ValueError
    

class Point:
    def __init__(self, x, y):
        if not (isarealnumber(x) and isarealnumber(y)): raise ValueError
        self.x = x
        self.y = y
        self.coordinates = (self.x, self.y)
        self.PO = self.dToB((0, 0))

    def __str__(self):
        return "({0}, {1})".format(self.x, self.y)

    def __repr__(self):
        return self.__str__()

    #  def __len__(self): return self.dToB((0, 0))
    #  unneeded since points do not have any property like length

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
        return self.info(return_Info=True, printInfo=False).replace("The polynom is", "A polynom with the equation", 1)

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
        return "This polynom is P(x) = " + str_editor(info)


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
        self.delta = (4*(self.delta_0 ** 3) - (self.delta_0 ** 2)) / (27*p*p)

    def __str__(self):
        return self.info().replace("polynom", "cubic", 1)

    def __repr__(self):
        return self.__str__()

    def the_roots(self):
        if self.delta_1 == 0 and self.delta_0 == 0:
            c_r = (-self.q / (3 * self.p))
            return (c_r, c_r, c_r)

p = Cubic(3,0,2,1)
print(p, type(p))
