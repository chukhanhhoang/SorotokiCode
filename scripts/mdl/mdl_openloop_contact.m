clr; cd;
%% 
L = 100;   % length of robot
% pcc
M = 20;     % number of modes/links
P = 4;      % points per link
% % chebyshev
% M = 8;     % number of modes/links
% P = 10;      %

N = M*P;  % number of discrete points on curve
H = 1/100; % timesteps
FPS = 30;  % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);

shp.E    = 0.05;     % Young's modulus in Mpa
shp.Nu   = 0.33;     % Poisson ratio
shp.Rho  = 1000e-12; % Density in kg/mm^3
shp.Zeta = 0.2;      % Damping coefficient

shp.Gvec = [0; 0; -150];

shp = shp.rebuild();

%% Object
xsph = 30;ysph=0; zsph = -30; rsph = 15;
sp=sSphere(xsph,ysph,-zsph,rsph);
sp2=sSphere(xsph,ysph,zsph,rsph-3);
%% Model with contact
mdl = cModel(shp,'Tstep',H,'Tsim',10);
mdl.planar = true;
mdl.object = sp;
mdl.pts_per_link = P;
mdl.object_center = [xsph; ysph; -zsph];
% mdl.constrained_points = [round(N/8)];
mdl.constraint_type = "pinned";

%%
% mdl.q0(1)    = 0.5;
% mdl.q0(2)    = -0.5;
% mdl.q0(3)    = -0.25;
mdl = mdl.simulate(); 
%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
[rig] = setupRig(M,L,Modes);
obj = Gmodel(sp2);
obj.Texture = diffuse(0.925);
obj.bake.render();

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end


%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
%    Y(:,ii) = pcc(X/L,ii,M); % pcc
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup rig
function [rig, gmdl] = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl');
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
     'SSSPower',0.005,'SSSRadius',5,'SSS',true);
 
gmdl = gmdl.bake.render(); 
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,base);
rig.g0 = SE3(roty(pi/2),zeros(3,1));

rig = rig.render();
end

