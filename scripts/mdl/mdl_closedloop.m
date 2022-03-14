clr; 
%% 
L = 100;  % length of robot
N = 40;   % number of discrete points on curve
M = 3;    % number of modes
H = 1/125; % timesteps
FPS = 30; % animation speed

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
shp.Zeta = 0.1;      % Damping coefficient

shp = shp.rebuild();
Theta_ = shp.get('ThetaEval');
Theta = pagemtimes(shp.Ba,Theta_);

%%
mdl = Model(shp,'Tstep',H,'Tsim',25);
mdl.gVec = [0;0;-9.81e3];
% mdl.gVec = [0;0;0];
% mdl = mdl.computeEL(mdl.q0);

%% find final config
gd = SE3(eye(3), [50;40;20]);

% qd = [-0.05;-0.1;-0.05;0.05;-0.02;0.01;0;0;0;0];
% p = shp.FK(qd);

%% controller
mdl.tau = @(M) Controller(M,gd);

%%
mdl.q0(1)   = 0.0;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
% figure;
% hold on;
[rig] = setupRig(M,L,Modes);

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    hold on;
%     plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
    scatter3(gd(1,4),gd(2,4),gd(3,4),100,col(2));
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end

%%
function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);
for ii = 1:M
   Y(:,ii) = legendre(X/L,ii-1); % chebyshev (basis function)
end
% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

function ret = intTheta(Theta,ds)
    ret = zeros(size(Theta));
    for i = 2:size(Theta,3)
        ret(:,:,i) = ret(:,:,i-1) + ds/2*(Theta(:,:,i)+Theta(:,:,i-1));
    end
end

%% setup controller
function [tau,error] = Controller(mdl,gd)
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    % 
    %tau        = zeros(n,1);
    g = SE3(mdl.Log.Phi(:,:,end),mdl.Log.p(:,end));
    
    J = mdl.Log.EL.J(:,:,end);
    
    k1 = 0.1;
    k2 = 0.1;
    lam1 = 2e-3;
    
    Kp = diag([k1,k1,k1,k2,k2,k2]);
    Xi = logmapSE3(g(:,:,end)\gd);
    Fu = Kp*tmapSE3(Xi)*isomse3(Xi);
    error = Fu; % just keep log
    dq = lam1*J.'*Fu;
    
    dV_dq = mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
    dVd_dq = -dq;
    tau = dV_dq-dVd_dq-2*mdl.Log.EL.M*mdl.Log.dq;
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