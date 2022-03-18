clr;

syms k1 k2 R d;

V = k1*d^2/2+ k2*(R-d)^2/2;

dVdd = simplify(diff(V,d))

d2Vd2d = diff(dVdd,d)

solve(dVdd==0,d)