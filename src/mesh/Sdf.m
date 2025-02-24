classdef Sdf
     
    properties
        sdf;
        BdBox;
        cmap = turbo;
        eps  = 1e-5;
        
        Node;      
        Element;
    end
%--------------------------------------------------------------------------    
    methods        
%---------------------------------------------------- Signed Distance Class     
        function obj = Sdf(fnc,varargin)
            obj.sdf = @(x) [fnc(x),fnc(x)];
            obj.cmap = turbo;
            for ii = 1:2:length(varargin)
                obj.(varargin{ii}) = varargin{ii+1};
            end
        end
%-------------------------------------------------------------------- union
        function r = plus(obj1,obj2)
            fnc = @(x) dUnion(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end
%--------------------------------------------------------------- difference        
        function r = minus(obj1,obj2)
            fnc = @(x) dDiff(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            r.BdBox = obj1.BdBox;
        end
%---------------------------------------------------------------- intersect        
        function r = mrdivide(obj1,obj2)
            fnc = @(x) dIntersect(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end
%----------------------------------------------------------- rotate X <-> Y              
        function r = transpose(obj1)
            fnc = @(x) obj1.sdf([x(:,2),x(:,1)]);
            r = Sdf(fnc);
            B = obj1.BdBox;
            r.BdBox = [B(3), B(4), B(1), B(2)];
        end
%----------------------------------------------------------- rotate X <-> Y              
        function r = rotate(obj1,varargin)
            if isempty(varargin), k = pi/2;
            else, k = varargin{1};
            end
            R = TwoDRotation(k);
            fnc = @(x) obj1.sdf((R*x.').');
            r = Sdf(fnc);
            B = obj1.BdBox;
            BB = [B(1),B(3);B(1),B(4);B(2),B(3);B(2),B(4)];
            
            r.BdBox = boxhull((R*BB.').').';
        end        
%----------------------------------------------------------- rotate X <-> Y              
        function r = smoothunion(obj1,obj2,varargin)
            if isempty(varargin), k = 1;
            else, k = varargin{1};
            end
                
            fnc = @(x) dSmoothUnion(obj1.sdf(x),obj2.sdf(x),k);
            r = Sdf(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end        
    end
%--------------------------------------------------------------------------       
methods (Access = public)
%--------------------------------------------------------- evalution of SDF
function d = eval(Sdf,x)
d = Sdf.sdf(x);
end
%------------------------------------------------- evalution of tangent SDF
function [T,N,B,Z] = normal(Sdf,x)
d = Sdf.sdf(x);
N = zeros(size(x,1),3);

if size(x,2) == 2
    n1 = (Sdf.sdf(x+repmat([Sdf.eps,0],size(x,1),1))-d)/Sdf.eps;
    n2 = (Sdf.sdf(x+repmat([0,Sdf.eps],size(x,1),1))-d)/Sdf.eps;
    %[Fy] = gradient(d(:,end));   
    %T = [Fy(:,1),Fy(:,2)];
    N(:,1) = n1(:,end);
    N(:,2) = n2(:,end);
else
    n1 = (Mesh.SDF(x+repmat([Sdf.eps,0,0],size(x,1),1))-d)/Sdf.eps;
    n2 = (Mesh.SDF(x+repmat([0,Sdf.eps,0],size(x,1),1))-d)/Sdf.eps;
    n3 = (Mesh.SDF(x+repmat([0,0,Sdf.eps],size(x,1),1))-d)/Sdf.eps;
    N = [n1(:,end),n2(:,end),n3(:,end)];
end

N = N./sqrt((sum((N.^2),2)));
T = ([0,1,0;-1,0,0;0,0,1]*N.').';
B = cross(T,N);
Z = atan2(N(:,2).*sign(-d(:,end)),N(:,1).*sign(-d(:,end)));
    
end
% %--------------------------------------------------------------------- show
% function [T, Z] = tangent(Sdf,x)
%  d = Sdf.sdf(x);
% N = zeros(size(x,1),3);
% 
% if size(x,2) == 2
%     n1 = (Sdf.sdf(x+repmat([Sdf.eps,0],size(x,1),1))-d)/Sdf.eps;
%     n2 = (Sdf.sdf(x+repmat([0,Sdf.eps],size(x,1),1))-d)/Sdf.eps;
%     %[Fy] = gradient(d(:,end));   
%     %T = [Fy(:,1),Fy(:,2)];
%     N(:,1) = n1(:,end);
%     N(:,2) = n2(:,end);
% else
%     n1 = (Mesh.SDF(x+repmat([Sdf.eps,0,0],size(x,1),1))-d)/Sdf.eps;
%     n2 = (Mesh.SDF(x+repmat([0,Sdf.eps,0],size(x,1),1))-d)/Sdf.eps;
%     n3 = (Mesh.SDF(x+repmat([0,0,Sdf.eps],size(x,1),1))-d)/Sdf.eps;
%     N = [n1(:,end),n2(:,end),n3(:,end)];
% end
% 
% N = N./sqrt((sum((N.^2),2)));
% T = ([0,1,0;-1,0,0;0,0,1]*N.').';
% Z = atan2(N(:,2).*sign(-d(:,end)),N(:,1).*sign(-d(:,end)));
%    
% end
%--------------------------------------------------------------------- show
function show(Sdf,Quality)
    
    if nargin < 2 
        if numel(Sdf.BdBox) < 6, Quality = 250;
        else, Quality = 50;
        end
    end
    
    x = linspace(Sdf.BdBox(1),Sdf.BdBox(2),Quality);
    y = linspace(Sdf.BdBox(3),Sdf.BdBox(4),Quality);
    
    if numel(Sdf.BdBox) < 6
        [X,Y] = meshgrid(x,y);

        D = Sdf.eval([X(:),Y(:)]);
        D = D(:,end);
        
        figure(101);
        surf(X,Y,reshape(D,[Quality Quality]),'linestyle','none');
        axis equal; hold on;
        
        if isempty(Sdf.Element)
            contour3(X,Y,reshape(D,[Quality Quality]),[0 0],'linewidth',...
                2.5,'Color','w');
        else
            patch('Vertices',Sdf.Node,'Faces',Sdf.Element,...
                'LineW',2.5,'EdgeColor','w');
        end
        
        colormap(Sdf.cmap);
        view(0,90);
    else
        z = linspace(Sdf.BdBox(5),Sdf.BdBox(6),Quality);
        [X,Y,Z] = meshgrid(x,y,z);

        D = Sdf.eval([X(:),Y(:),Z(:)]);
        D = D(:,end);
        D(D>0) = NaN;
        figure(101);
        C = reshape(D,[Quality Quality,Quality]);
        scatter3(X(:),Y(:),Z(:),85,C(:),'Marker','.');
        axis equal; hold on;

        colormap(Sdf.cmap);
    end
end
%------------------------------------------------------------- show contour
function showcontour(Sdf,Quality)
    if nargin < 2 
        if numel(Sdf.BdBox) < 6
            Quality = 250;
        else
            Quality = 25;
        end
    end
    
    x = linspace(Sdf.BdBox(1),Sdf.BdBox(2),Quality);
    y = linspace(Sdf.BdBox(3),Sdf.BdBox(4),Quality);
    
    if numel(Sdf.BdBox) < 6
        [X,Y] = meshgrid(x,y);

        D = Sdf.eval([X(:),Y(:)]);
        D = D(:,end);
        
        %figure(101);
        contour3(X,Y,reshape(D,[Quality Quality]),[0 0],'linewidth',...
            2,'Color','k'); hold on;
        
        D(D<=0) = 0;
        D(D>1e-6)  = NaN;
        
        cplane(X,Y,reshape(D,[Quality Quality]));
        %axis equal; 
        
        
        colormap(Sdf.cmap);
        %view(0,90);
        %axis off;
        caxis([-1 1]);
    end
end
end
    
methods (Access = private)

end
end

function [R] = TwoDRotation(x)
R = [cos(x),sin(x);-sin(x),cos(x)];
end
