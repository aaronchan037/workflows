import gmpy2
from gmpy2 import is_prime, isqrt, powmod, invert

def monty_gen_basept(m, A):
    '''
    蒙哥马利曲线生成基点
    '''
    # 定义有限域
    p = 2**m - 1
    while not is_prime(p):
        p -= 2

    # 转换A到有限域
    A = int(A) % p

    # 查找生成元
    x = 1
    while True:
        y_squared = (pow(x, 3, p) + A * pow(x, 2, p) + x) % p
        if powmod(y_squared, (p - 1) // 2, p) == 1:  # 判断是否是平方
            y = isqrt(y_squared)
            G = (x, y)
            break
        x += 1

    # 查找基点
    x = 1
    while True:
        y_squared = (pow(x, 3, p) + A * pow(x, 2, p) + x) % p
        if powmod(y_squared, (p - 1) // 2, p) == 1:  # 判断是否是平方
            y = isqrt(y_squared)
            P = (x, y)
            if not is_point_on_curve(P, A, p) or not has_order(P, A, p):
                continue
            break
        x += 1

    print(f"Monty Generator: {G}")
    print(f"Monty Base Point: {P}")

def is_point_on_curve(P, A, p):
    x, y = P
    return pow(y, 2, p) == (pow(x, 3, p) + A * pow(x, 2, p) + x) % p

def has_order(P, A, p):
    n = 8 * (p + 1 - (2 * (p + 1) // 8))
    r = n // 8
    if not is_prime(r):
        return False
    x, y = P
    P8 = point_multiply(P, 8, A, p)
    Pr = point_multiply(P, r, A, p)
    P2r = point_multiply(P, 2*r, A, p)
    P4r = point_multiply(P, 4*r, A, p)
    Pn = point_multiply(P, n, A, p)
    return P8 != (0, 0) and Pr != (0, 0) and P2r != (0, 0) and P4r != (0, 0) and Pn == (0, 0)

def point_double(P, A, p):
    x1, y1 = P
    if y1 == 0:
        return (0, 0)
    lambda_ = (3 * pow(x1, 2, p) + 2 * A * x1 + 1) * invert(2 * y1, p) % p
    x3 = (pow(lambda_, 2, p) - 2 * x1) % p
    y3 = (lambda_ * (x1 - x3) - y1) % p
    return (x3, y3)

def point_add(P1, P2, A, p):
    x1, y1 = P1
    x2, y2 = P2
    if P1 == (0, 0):
        return P2
    if P2 == (0, 0):
        return P1
    if x1 == x2 and y1 != y2:
        return (0, 0)
    if P1 == P2:
        return point_double(P1, A, p)
    lambda_ = (y2 - y1) * invert(x2 - x1, p) % p
    x3 = (pow(lambda_, 2, p) - x1 - x2) % p
    y3 = (lambda_ * (x1 - x3) - y1) % p
    return (x3, y3)

def point_multiply(P, n, A, p):
    result = (0, 0)
    addend = P
    while n > 0:
        if n % 2 == 1:
            result = point_add(result, addend, A, p)
        addend = point_double(addend, A, p)
        n //= 2
    return result
