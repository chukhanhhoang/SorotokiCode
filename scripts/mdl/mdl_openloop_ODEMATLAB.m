clr; cd;
%% 
L = 100;   % length of robot
M = 4;     % number of modes
N = M*10;  % number of discrete points on curve
H = 1/75; % timesteps
FPS = 30;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature
%%
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);

shp.E    = 0.05;     % Young's modulus in Mpa
shp.Nu   = 0.33;     % Poisson ratio
shp.Rho  = 1000e-12; % Density in kg/mm^3
shp.Zeta = 0.1;      % Damping coefficient

shp.Gvec = [0; 0; -9.81];

shp = shp.rebuild();

%%
mdl = Model(shp,'Tstep',H,'Tsim',15);
% get evals
mdl.Theta = mdl.Shapes.get('ThetaEval');
mdl.Xi0 = mdl.Shapes.get('Xi0Eval');

%%
mdl.q0(1)    = 0;

tic
[t,x] = ode45(@(t,x) dynamics(t,x,mdl),[0 mdl.Tsim],[mdl.q0;zeros(size(mdl.q0))]);
toc
% mdl = mdl.simulate(); 
%% 

mdl.Log.t = t;
mdl.Log.q = x(:,1:size(x,2)/2);

figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
[rig] = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end

function dxdt = dynamics(t,x,Model)
% x = [q;dq] state of the robot
    Q = x(1:(numel(x)/2));
    dQ = x((numel(x)/2+1):end);

    [M_,C_,K_,R_,G_,...
                p_,Phi_,J_,Vg_,Kin_] = computeLagrangianFast_mex(...
                Q,dQ,... 
                Model.Shapes.ds,...   
                Model.p0,... 
                Model.Phi0,...
                Model.Xi0,... 
                Model.Theta,...
                Model.Shapes.Ba,... 
                Model.Shapes.Ktt,...
                Model.Shapes.Mtt,...     
                Model.Shapes.Zeta,...
                Model.Shapes.Gvec);
    Minv = M_\eye(numel(Q));
    dxdt = [dQ;Minv*(- C_*dQ - K_*Q - R_*dQ - G_)];
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
rig.g0 = SE3(roty(-pi),zeros(3,1));

rig = rig.render();
end