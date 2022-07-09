# getting Edwards form

from math import gcd
from Crypto.Util.number import *


########## WORKS IF THE NUMBER OF POINTS OF THE CURVE IS A MULTIPLE OF 4
########## Otherwise Montgomery and Edwards form do not exist


def ison(C, P):
    """
    Verification points are on the curve
    """
    c, d, p = C
    u, v = P
    return (u**2 + v**2 - cc * (1 + d * u**2*v**2)) % p == 0

def a_and_b(u1,u2,v1,v2):
    """
    Helper function used to simplify calculations
    """
    a12 = u1**2 - u2**2 + v1**2 - v2**2
    b12 = u1**2 * v1**2 - u2**2 * v2**2
    return a12, b12

def find_modulus(u1,u2,u3,u4,v1,v2,v3,v4):
    """
    Compute the modulus from four points
    """
    a12, b12 = a_and_b(u1,u2,v1,v2)
    a13, b13 = a_and_b(u1,u3,v1,v3)
    a23, b23 = a_and_b(u2,u3,v2,v3)
    a24, b24 = a_and_b(u2,u4,v2,v4)

    p_almost = gcd(a12*b13 - a13*b12, a23*b24 - a24*b23)

    for i in range(2,1000):
        if p_almost % i == 0:
            p_almost = p_almost // i

    return p_almost

def c_sq_d(u1,u2,v1,v2,p):
    """
    Helper function to computer c^2 d
    """
    a1,b1 = a_and_b(u1,u2,v1,v2)
    return a1 * pow(b1,-1,p) % p

def c(u1,u2,v1,v2,p):
    """
    Compute c^2, d from two points and known modulus
    """
    ccd = c_sq_d(u1,u2,v1,v2,p)
    cc = (u1**2 + v1**2 - ccd*u1**2*v1**2) % p
    d = ccd * pow(cc, -1, p) % p
    return cc, d


P = (57965690875794075971012488152230149163602256508707811594337978733750915468806, 16746479600411735837174391941067872880449998685866235929656201743399769433149)
Q = (69936234690532721081726160313018549398162634031107454554647599731204812649342, 60370347819144291939385104514502372001122251366410414587034042305240897452453)
sP = (68863793597568829791091117378551243306095213130450629884577755357985078854249, 57338665993204791888869144252570837621696803230095730917988989291398052731921)
tQ = (21561479239953076675545350206846696529573997309059024341216968862074706660695, 45298495518669860523224015532327258014763911346364900700009286850163858621557)

u1, v1 = P
u2, v2 = Q
u3, v3 = sP
u4, v4 = tQ

p = find_modulus(u1,u2,u3,u4,v1,v2,v3,v4)
cc, d = c(u1,u2,v1,v2,p)

C = cc, d, p
assert ison(C, P)
assert ison(C, Q)
assert ison(C, sP)
assert ison(C, tQ)

print(f'Found curve parameters')
print(f'p = {p}')
print(f'c^2 = {cc}')
print(f'd = {d}')

# Found curve
# p = 903968861315877429495243431349919213155709
# cc = 495368774702871559312404847312353912297284
# d = 540431316779988345188678880301417602675534


F = GF(p)

c = F(cc).sqrt()
x1, y1 = P
x2, y2 = Q
x3, y3 = sP
x4, y4 = tQ

R.<x,y> = PolynomialRing(F)
g = (x^2 + y^2 - cc * (1 + d * x^2*y^2))

# Check the mapping worked!
assert g(x=x1, y=y1) == 0
assert g(x=x2, y=y2) == 0
assert g(x=x3, y=y3) == 0
assert g(x=x4, y=y4) == 0
###### ------------------------------------------------------------------ From Edwards to Montgomery
# Scale: x,y,d to remove c:
# x^2 + y^2 = c^2 * (1 + d * x^2*y^2)
# to:
# x^2 + y^2 = (1 + d * x^2*y^2)

d = F(d) * F(cc)^2
x1, y1 = F(x1) / F(c),  F(y1) / F(c)
x2, y2 = F(x2) / F(c),  F(y2) / F(c)
x3, y3 = F(x3) / F(c),  F(y3) / F(c)
x4, y4 = F(x4) / F(c),  F(y4) / F(c)

h = (x^2 + y^2 - (1 + d * x^2*y^2))

# Check the mapping worked!
assert h(x=x1, y=y1) == 0
assert h(x=x2, y=y2) == 0
assert h(x=x3, y=y3) == 0
assert h(x=x4, y=y4) == 0

# Convert from Edwards to Mont. 
# https://safecurves.cr.yp.to/equation.html
def ed_to_mont(x,y):
    u = F(1 + y) / F(1 - y)
    v = 2*F(1 + y) / F(x*(1 - y))
    return u,v

#######################################------------------- getting Weierstrass from Montgomery


u1, v1 = ed_to_mont(x1, y1)
u2, v2 = ed_to_mont(x2, y2)
u3, v3 = ed_to_mont(x3, y3)
u4, v4 = ed_to_mont(x4, y4)

e_curve = 1 - F(d)
A = (4/e_curve - 2)
B = (1/e_curve)

# Mont. curve: Bv^2 = u^3 + Au^2 + u
R.<u,v> = PolynomialRing(ZZ)
f = B*v^2 - u^3 - A* u^2 - u

# Check the mapping worked!
assert f(u=u1, v=v1) == 0
assert f(u=u2, v=v2) == 0
assert f(u=u3, v=v3) == 0
assert f(u=u4, v=v4) == 0

# Convert from Mont. to Weierstrass
# https://en.wikipedia.org/wiki/Montgomery_curve
a = F(3 - A^2) / F(3*B^2)
b = (2*A^3 - 9*A) / F(27*B^3)
E = EllipticCurve(F, [a,b])

# https://en.wikipedia.org/wiki/Montgomery_curve
def mont_to_wei(u,v):
    t = (F(u) / F(B)) + (F(A) / F(3*B))
    s = (F(v) / F(B))
    return t,s

X1, Y1 = mont_to_wei(u1, v1)
X2, Y2 = mont_to_wei(u2, v2)
X3, Y3 = mont_to_wei(u3, v3)
X4, Y4 = mont_to_wei(u4, v4)

P = E(X1, Y1)
Q = E(X2, Y2)
sP = E(X3, Y3)
tQ = E(X4, Y4)


print(a)
print(b)