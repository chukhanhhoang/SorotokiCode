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
xs = 30; ys = 7; zs = -45; rs = 15;
sphere_pos = [xs;ys;zs];

%% controller
mdl.tau = @(M) Controller(M,sphere_pos,rs);

%%
mdl.q0(1)   = 0;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
[rig] = setupRig(M,L,Modes);

[X,Y,Z] = sphere();
X = rs*X+xs;
Y = rs*Y+ys;
Z = rs*Z+zs;

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)

    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    hold on;
    surf(X,Y,Z);hold on;
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(0,20)
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
function tau = Controller(mdl,sphere_pos,rs)
n = size(mdl.Log.q,1);
t = mdl.Log.t;

% % Sphere position, radius
% xs = 30; ys = 0; zs = -20; rs = 10;
% sphere_pos = [xs;ys;zs];
stiffness = 1e-3;

% Init
tau        = zeros(n,1);

[g,J] = mdl.Shapes.string(mdl.Log.q);

position = reshape(g(1:3,4,:),3,[]);

for i = 1:size(position,2)
    vector_from_sphere = position(:,i)-sphere_pos;
    if norm(vector_from_sphere) <= rs
        body_force = [zeros(3,1);g(1:3,1:3,i).'*stiffness*(1-1/norm(vector_from_sphere))*vector_from_sphere];
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

function plotAxes(hAx)
    axis(hAx,'equal')
    %Get X, Y and Z data for plotting the axes...
    X_range = hAx.XLim(2) - hAx.XLim(1);
    X_start = hAx.XLim(1);
    X_delta = X_range/20;
    Y_delta = (hAx.YLim(2) - hAx.YLim(1))/20;
    Y_start = hAx.YLim(1);
    Z_delta = (hAx.ZLim(2) - hAx.ZLim(1))/20;
    Z_start = hAx.ZLim(1);
    X_Line = line(hAx,[X_start+X_delta X_start+X_delta*5],[Y_start+Y_delta Y_start+Y_delta],[Z_start+Z_delta Z_start+Z_delta]); % x Line
    Y_Line = line(hAx,[X_start+X_delta X_start+X_delta],[Y_start+Y_delta Y_start+Y_delta*5],[Z_start+Z_delta Z_start+Z_delta]); % Y Line
    Z_Line = line(hAx,[X_start+X_delta X_start+X_delta],[Y_start+Y_delta Y_start+Y_delta],[Z_start+Z_delta Z_start+Z_delta*5]); %Z Line
    X_text = text(hAx,X_start+X_delta*6,Y_start+Y_delta,Z_start+Z_delta,'x');
    Y_text = text(hAx,X_start+X_delta,Y_start+Y_delta*6,Z_start+Z_delta,'y');
    Z_text = text(hAx,X_start+X_delta,Y_start+Y_delta,Z_start+Z_delta*6,'z');
end