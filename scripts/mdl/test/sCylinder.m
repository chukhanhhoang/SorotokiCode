function sdf = sCylinder(xc,yc,zc,r,l)
if nargin < 2
    r = xc; 
    xc = 0;
    yc = 0;
    zc = 0;
    l = 1;
end

sdf = Sdf(@(P) sdfCylinder(P,xc,yc,zc,r,l));

% % r = 2*r;

sdf.BdBox = [xc-r-1e-6,xc+r+1e-6,....
             yc-r-1e-6,yc+r+1e-6,...
             zc-r-1e-6,zc+r+1e-6];
         
[sdf.Node,sdf.Element] = generateNodeSet(xc,yc,zc,r,l,30);
         
end

function d = sdfCylinder(P,xc,yc,zc,r,l)
n = size(P,1);
d = zeros(n,1);
for i = 1:n
    if P(3,i)<= zc+l/2 && P(3,i)>= zc-l/2
        d(i) = min(sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
    end
end
% d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2+(P(:,3)-zc).^2)-r;
d = [d,d];
end

function [V,F] = generateNodeSet(xc,yc,zc,r,l,N)
th = linspace(-pi,pi,N).';

x1 = r*cos(th) + xc;
y1 = r*sin(th) + yc;
z1 = x1*0;

x2 = x1;
y2 = y1;
z2 = r*cos(th) + zc;

x3 = x1*0;
y3 = r*sin(th) + xc;
z3 = r*cos(th) + zc;

V = [x1,y1,z1;x2,y2,z2;x3,y3,z3];
S = [(1:N-1).',(2:N).'];
F = [S;S+N;S+2*N];
%F = [1:N;N+1:2*N;2*N+1:3*N];
end