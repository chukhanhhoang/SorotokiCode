clr;

syms gQ [7 1];
syms xi [6 1];
syms lambda [7 1];
syms xi0 [6 1];
syms K [6 6];
syms r_o [3 1]

A = [0     -xi(1)   -xi(2)   -xi(3);
    xi(1)   0        xi(3)   -xi(2);
    xi(2)  -xi(3)    0        xi(1);
    xi(3)   xi(2)   -xi(1)    0];

f = [0.5*A*gQ(1:4);quat2rot(gQ(1:4))*xi(4:6)];

H = lambda.'*f - (xi-xi0).' * K * (xi-xi0);

hatH = H - sumsqr(gQ(5:7)-r_o);

dHhat_dgQ = matlabFunction(partialdiff(hatH,gQ),'Vars',{gQ,xi,lambda,xi0,r_o,K});
dHhat_dxi = matlabFunction(partialdiff(hatH,xi),'Vars',{gQ,xi,lambda,xi0,r_o,K});

function ret = sumsqr(x)
    ret = x.'*x;
end

function ret = partialdiff(f,x)
    n = length(x);
    ret = sym('temp',[1,n]);
    for i = 1:n
        ret(i)=diff(f,x(i));
    end
end

function jacobian = calJacobian(f,x)
    jacobian = sym('temp',[length(f),length(f)]);

    for i = 1:length(f)
        jacobian(i,:) = partialdiff(f(i),x);
    end
end