clr; 
%% 
L = 100;  % length of robot
N = 50;   % number of discrete points on curve
M = 5;    % number of modes
H = 1/125; % timesteps
FPS = 125; % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);

shp.E    = 2.00;     % Young's modulus in Mpa
shp.Nu   = 0.49;     % Poisson ratio
shp.Rho  = 1000e-12; % Density in kg/mm^3
shp.Zeta = 0.15;      % Damping coefficient

shp = shp.rebuild();

%% Init model
mdl = Model(shp,'Tstep',H,'Tsim',10);
mdl.gVec = [0;0;-9.81e3];

% Sphere position, radius
xs = 10; ys = 7; zs = -65; rs = 20;
sphere = sSphere(xs+rs,ys,zs+rs,rs); % sphere
cube = sCube(xs-rs,xs+rs,ys-rs,ys+rs,zs-rs,zs+rs); % cube centered at xs,ys,zs
sdf = sphere + cube;
msh = Mesh(sdf,'NElem',1e3);
sphere_gmodel = Gmodel(sdf);

%% controller
mdl.tau = @(M) Controller(M,msh);

%%
mdl.q0(1)   = 0;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
sphere_gmodel.bake()
[rig] = setupRig(M,L,Modes);


for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    hold on;
    sphere_gmodel.render();
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(45,-15)
    drawnow();
end


%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

%% setup controller
function tau = Controller(mdl,obj_mesh)
n = size(mdl.Log.q,1);
t = mdl.Log.t;

% % Sphere position, radius
% xs = 30; ys = 0; zs = -20; rs = 10;
% sphere_pos = [xs;ys;zs];
stiffness = 5e-3;
damping_base = 2e-5;
% Init
tau        = zeros(n,1);

[g,J] = mdl.Shapes.string(mdl.Log.q);


position = reshape(g(1:3,4,:),3,[]);
[N,dist] = obj_mesh.computeN(position.');

for i = 1:size(position,2)
    if dist(i) <= 0
        body_velo_twist = J(:,:,i)*mdl.Log.dq;
        spatial_velo_twist = Admap(g(1:3,1:3,i),g(1:3,4,i))*body_velo_twist;
        spatial_velo = isomse3(spatial_velo_twist)*[g(1:3,4,i);0];
        dd = N(i,:)*spatial_velo(1:3);
        spatial_damping_force = -damping_base*((-dist(i))^1.1)*dd*N(i,:).';
        spatial_stiffness_force=(-dist(i))*stiffness*N(i,:).';
        
        body_force = [zeros(3,1);g(1:3,1:3,i).'*(spatial_stiffness_force+spatial_damping_force)];
        tau = tau + J(:,:,i).' * body_force;
    end
end
% tau(1)     = 9*smoothstep(t)*sin(3*t);
% tau(n/2+1) = 9*smoothstep(t)*cos(t);
end

%% setup rig
function [rig, gmdl] = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl');
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
     'SSSPower',0.005,'SSSRadius',5,'SSS',true);
 
N = 200;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,mateplastic);
rig.g0 = SE3(roty(-pi),zeros(3,1));

rig = rig.render();
end