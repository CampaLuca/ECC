def tonelli_shanks(n, p):
    if pow(n, int((p-1)//2), p) == 1:
            s = 1
            q = int((p-1)//2)
            while True:
                if q % 2 == 0:
                    q = q // 2
                    s += 1
                else:
                    break
            if s == 1:
                r1 = pow(n, int((p+1)//4), p)
                r2 = p - r1
                return r1, r2
            else:
                z = 2
                while True:
                    if pow(z, int((p-1)//2), p) == p - 1:
                        c = pow(z, q, p)
                        break
                    else:
                        z += 1
                r = pow(n, int((q+1)//2), p)
                t = pow(n, q, p)
                m = s
                while True:
                    if t == 1:
                        r1 = r
                        r2 = p - r1
                        return r1, r2
                    else:
                        i = 1
                        while True:
                            if pow(t, 2**i, p) == 1:
                                break
                            else:
                                i += 1
                        b = pow(c, 2**(m-i-1), p)
                        r = r * b % p
                        t = t * b ** 2 % p
                        c = b ** 2 % p
                        m = i
    else:
        return False