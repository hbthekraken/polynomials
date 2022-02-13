from polynomials_extended_011 import Line

digits = '0123456789'

warn_message = str("Something went wrong, check your equations please. There might be:" +
                   "\n-equal signs \'=\' more than one,\n-unexpected terms such as \'xy\'," +
                   "\n-or some unexpected/wrong symbols/notations being used!")


def find_x_y(alpha, a, b, beta, c, d):
    l1 = Line((a, b, -alpha))
    l2 = Line((c, d, -beta))
    return l1.interception_with(l2)
    """
    if l1.wrt_l2(l2) == 'C':
        return None
    elif l1.wrt_l2(l2) == 'P':
        return (None, None)
    else:
        if a != 0:
            y = (beta - ((c * alpha) / a)) / (d - ((c * b) / a))
            x = (alpha - (b * y)) / a
            return (x, y)
        else:
            y = alpha / b
            x = (beta - d * y) / c
            return (x, y)
    """


def isthereanother(str1, charss):
    if len(charss) == 1:
        return True
    if len(str1) == 1:
        return True
    for c in charss:
        str1 = str1.replace(c, "")
    return len(str1) != 0


def check_the_eq(_str):
    if _str.count("=") != 1:
        return False
    if "xy" in _str or "yx" in _str:
        return False
    if isthereanother(_str, digits + "xy .+-="):
        return False
    return True


def coefficient_extractor(l):
    coef = {
        "x": 0,
        "y": 0,
        "c": 0
    }
    for i in l:
        if "x" in i:
            if "-x" in i:
                coef.update({"x": coef["x"] - 1.0})
            elif i.strip("x") != "":
                coef.update({"x": coef["x"] + float(i.strip("x"))})
            else:
                coef.update({"x": coef["x"] + 1.0})
        elif "y" in i:
            if "-y" in i:
                coef.update({"y": coef["y"] - 1.0})
            elif i.strip("y") != "":
                coef.update({"y": coef["y"] + float(i.strip("y"))})
            else:
                coef.update({"y": coef["y"] + 1.0})

        elif len(i) == 0:
            continue
        else:
            #  print("i = ->{0}<-, {1}".format(i, type(i)))
            coef.update({"c": coef["c"] + float(i)})

    return coef


def equation_rearranger(s):
    s = str(s)
    ioes = s.find("=")
    l_of_eq = s[0:ioes]
    r_of_eq = s[ioes + 1:]
    l_coefs = coefficient_extractor(l_of_eq.split())
    r_coefs = coefficient_extractor(r_of_eq.split())
    total_coefs = {
        "x": l_coefs["x"] - r_coefs["x"],
        "y": l_coefs["y"] - r_coefs["y"],
        "c": l_coefs["c"] - r_coefs["c"]
    }
    return total_coefs


#  assuming that the variabes will be given as x and y, in lowercase
def linear_equation_solver(str1, str2):
    if not (isinstance(str1, str) and isinstance(str2, str)):
        raise ValueError
    if not (check_the_eq(str1) and check_the_eq(str2)):
        raise ValueError(warn_message)
    str1 = str1.replace(" ", "").replace("+", " ").replace("-", " -")
    t_c1 = equation_rearranger(str1)
    str2 = str2.replace(" ", "").replace("+", " ").replace("-", " -")
    t_c2 = equation_rearranger(str2)
    return find_x_y(-t_c1["c"], t_c1["x"], t_c1["y"], -t_c2["c"], t_c2["x"], t_c2["y"])


if __name__ == "__main__":
    print("Please enter equations that comprise only "
          + "x, y and constant terms all of whose coefficients are real numbers.")
    t = True
    while t:
        eq1 = input("Equation one: ")
        eq2 = input("Equation two: ")

        a = linear_equation_solver(eq1, eq2)
        print("The solution for this equation system, (x, y), is:", a)
        temp_input = input("Do you want to try again? Then type y and press enter, otherwise just press enter.\n")
        if temp_input.lower().strip() == "y":
            continue
        else:
            t = False
