"""
Credits to: https://gist.github.com/pqlx/
"""

from sage.all import *

def check_pohlig_hellman(curve, generator=None):
    """
    The Pohlig-Hellman algorithm allows for quick (EC)DLP solving if the order of the curve is smooth,
    i.e its order is a product of multiple (small) primes.
    The best general purpose algorithm for finding a discrete logarithm is the Baby-step giant-step
    algorithm, with a running time of O(sqrt(n)).
    If the order of the curve (over a finite field) is smooth, we can however solve the (EC)DLP
    algorithm by solving the (EC)DLP for all the prime powers that make up the order, then using the
    Chinese remainder theorem to compute the (EC)DLP solution to our original order.
    If no generator is provided, we assume a cofactor of 1, and thus a generator subgroup order equal to the curve order.
    """


    if generator:
        order = generator.order()
    else:
        order = curve.order()
    factorization = factor(order)

    # baby-step giant-step complexity is O(sqrt(n))
    naive_complexity = order.nth_root(2, True)[0] + 1

    # for an order N = (p_0)^(e_0) * (p_1)^(e_1) * ... * (p_k)^(e^k) with p prime and e natural the complexity is:
    # O( \sum_{i=0}^k e_i ( log_2(n) + sqrt((p_i)) ) )

    logn = log(order, 2)

    pohlig_hellman_complexity  = sum( y * (logn + x.nth_root(2, True)[0]) for x, y in factorization)

    return (pohlig_hellman_complexity, naive_complexity)


def check_smart(curve):
    """
    The Smart attack allows for solving the ECDLP in linear time, given that the trace of Frobenius of the curve equals one,
    i.e the order of the curve is equal to the order of the finite field over which the elliptic curve is defined.
    If this is the case, we can create a group isomorphism from our curve E over Fp to the finite field Fp (which preserves addition).
    Solving the discrete 'log' problem in a finite field is easy, when our group is an additive one.
    Finding `m` such that mQ = P given P, Q is hard in an elliptic curve defined over a finite field, but easy over a finite field an sich:
    m = P * modinv(Q)
    """

    return curve.trace_of_frobenius() == 1

def check_mov(curve):
    """
    The MOV attack (Menezes, Okamoto, Vanstone) allows for solving the ECDLP in subexponential time, given that the trace of frobenius of the curve equals zero,
    i.e the curve is supersingular.
    If this is the case we can reduce the problem to the discrete logarithm problem in the finite field over which the curve is defined. Subexponential attacks are known.
    This differs from the smart attack in the sense that we have to solve the actual multiplicative discrete log problem Q^m = P,
    instead of the additive discrete log problem mQ = P
    """

    return curve.trace_of_frobenius() == 0

def check(curve, generator=None):

    def print_red(p):
        print("\u001b[41m" + str(p) + "\u001b[0m")

    ph = check_pohlig_hellman(curve, generator)

    if ph[0] < ph[1]:
        quotient = round( float(ph[1]/ph[0]), 2) - 1

        logs = [round(float(log(x, 2)), 2) for x in ph]

        print_red(f"Pohlig-Hellman can make ECDLP solving {quotient} times faster!")
        print(f"O(2^{logs[1]}) -> O(2^{logs[0]})")

    else:
        print("Pohlig-Hellman cannot improve ECDLP solving speed.")

    smart = check_smart(curve)

    if smart:
        print_red("Smart's attack can solve ECDLP in linear time!")
    else:
        print("Smart's attack does not apply.")

    mov = check_mov(curve)

    if mov:
        print_red("MOV attack can solve ECDLP in subexponential time!")
    else:
        print("MOV attack does not apply.")

"""
example usage:
p = 310717010502520989590157367261876774703
a = 2
b = 3
p = 0xa15c4fb663a578d8b2496d3151a946119ee42695e18e13e90600192b1d0abdbb6f787f90c8d102ff88e284dd4526f5f6b6c980bf88f1d0490714b67e8a2a2b77
a = 0x5e009506fcc7eff573bc960d88638fe25e76a9b6c7caeea072a27dcd1fa46abb15b7b6210cf90caba982893ee2779669bac06e267013486b22ff3e24abae2d42
b = 0x2ce7d1ca4493b0977f088f6d30d9241f8048fdea112cc385b793bce953998caae680864a7d3aa437ea3ffd1441ca3fb352b0b710bb3f053e980e503be9a7fece
c = EllipticCurve(GF(p), [a, b])
check(c)
"""
