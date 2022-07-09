#Q = k*G
# We can solve the discrete logarithm by using the following algorithm if 

'''
Most cryptosystems based on elliptic curves can be broken if you can solve the discrete logarithm problem, that is, given the point P
and rP, find the integer r.

The MOV attack uses a bilinear pairing, which (roughly speaking) is a function e
that maps two points in an elliptic curve E(Fq) to a element in the finite field Fqk, where k is the embedding degree associated with the curve. The bilinearity means that e(rP,sQ)=e(P,Q)rs for points P,Q. Therefore, if you want to compute the discrete logarithm of rP, you can instead compute u=e(P,Q) and v=e(rP,Q) for any Q. Due to bilinearity, we have that v=e(P,Q)r=ur. Now you can solve the discrete logarithm in Fqk (given ur and u, find r) in order to solve the discrete logarithm in the elliptic curve!

Usually, the embedding degree k
is very large (the same size as q), therefore transfering the discrete logarithm to Fqk won't help you. But for some curves the embedding degree is small enough (specially supersingular curves, where k<=6), and this enables the MOV attack. For example, a curve with a 256-bit q usually offers 128 bits of security (i.e. can be attacked using 2128 steps); but if it has an embedding degree 2, then we can map the discrete logarithm to the field Fq2
which offers only 60 bits of security.
In practice the attack can be simply avoided by not using curves with small embedding degree; standardized curves are safe. Since pairings also have many constructive applications, it is possible to carefully choose curves where the cost of attacking the elliptic curve itself or the mapped finite field is the same.'''


def movAttack(G, Q):
    k = 1
    while (p**k - 1) % E.order():
        k += 1

    Ee = EllipticCurve(GF(p**k, 'y'), [a, b])

    R = Ee.random_point()
    m = R.order()
    d = gcd(m, G.order())
    B = (m // d) * R

    assert G.order() / B.order() in ZZ
    assert G.order() == B.order()

    Ge = Ee(G)
    Qe = Ee(Q)

    n = G.order()
    alpha = Ge.weil_pairing(B, n)
    beta = Qe.weil_pairing(B, n)

    print('Computing log...')
    nQ = beta.log(alpha)
    return nQ