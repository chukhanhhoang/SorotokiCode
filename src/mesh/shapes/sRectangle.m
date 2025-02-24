function sdf = sRectangle(x1,x2,y1,y2)
if nargin < 2
    if numel(x1) == 1
        a = x1;
        x1 = -a;
        x2 = a; 
        y1 = -a;
        y2 = a;
    else
        X = x1;
        x1 = X(1);
        x2 = X(2);
        y1 = X(3);
        y2 = X(4);
    end
end

eps = 1e-2*norm([x1;x2;y1;y2]);

sdf = Sdf(@(P) sdfRectangle(P,x1,x2,y1,y2));
sdf.BdBox = [x1-eps,x2+eps,y1-eps,y2+eps];
sdf.Node = [x1,y1;x2,y1;x2,y2;x1,y2;x1,y1];
sdf.Element = [(1:4).',(2:5).'];

end

function d = sdfRectangle(P,x1,x2,y1,y2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
d = [d,max(d,[],2)];
end