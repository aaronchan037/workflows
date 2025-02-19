import gmpy2
from gmpy2 import mpz, invert, powmod

def is_point_on_curve(P, A, p):
    x, y = P
    return powmod(y, 2, p) == (powmod(x, 3, p) + A * powmod(x, 2, p) + x) % p

def point_double(P, A, p):
    x1, y1 = P
    if y1 == 0:
        return (0, 0)
    lambda_ = (3 * powmod(x1, 2, p) + 2 * A * x1 + 1) * invert(2 * y1, p) % p
    x3 = (powmod(lambda_, 2, p) - 2 * x1) % p
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
    x3 = (powmod(lambda_, 2, p) - x1 - x2) % p
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

def verify_point(P, A, p, r):
    # 检查点是否在曲线上
    if not is_point_on_curve(P, A, p):
        print("点不在曲线上")
        return False

    # 检查点的阶数是否为 r
    O = point_multiply(P, r, A, p)
    if O != (0, 0):
        print("点的阶数不正确")
        return False

    print("点验证通过！")
    return True

# 示例参数
p = mpz(2**255 - 19)  # Curve25519 的素数
A = mpz(486662)
# 使用 Curve25519 的标准基点
P = (mpz(9), mpz(1478161944758954479102059356840998688726460613461607578403730737476650733120107129202958574023890561561637560978246667051))
r = mpz(2**252 + 27742317777372353535851937790883648493)  # Curve25519 的阶数

# 验证点
verify_point(P, A, p, r)
