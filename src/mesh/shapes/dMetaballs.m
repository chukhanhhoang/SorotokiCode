function d = dMetaballs(P,xc,yc,R,w)
d = zeros(length(P),1);

for ii = 1:length(xc)
    r = sqrt((P(:,1)-xc(ii)).^2+(P(:,2)-yc(ii)).^2);
    d = PiecewiseGaussian(r,w,R) + d; 
end
d = d+1e-3;
d = [d,d];
end

%-------------------------------------------------------------------------%
function d = PiecewiseGaussian(r,a,b)
 I1 = (r<=b/3) & (r>=0);
 I2 = (r<=b) & (r>b/3);
 I3 = (r>b);
 
 d = zeros(length(r),1);
 d(I1) = -a*(1-3*(b.^-2).*r(I1).^2);
 d(I2) = -(3*a/2).*(1-r(I2)/b).^2;
 d(I3) = eps;
end
