clr; 
%% 
L = 100;  % length of robot
N = 30;   % number of discrete points on curve
M = 5;    % number of modes
H = 1/125; % timesteps
FPS = 25; % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature
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
mdl = Model(shp,'Tstep',H,'Tsim',3);
mdl.gVec = [0;0;-4.81e3];

% Sphere position, radius
xs = 20; ys = 0; zs = -45; rs = 20;
sphere = sSphere(xs,ys,zs,rs); % sphere
sphere_gmodel = Gmodel(sphere);

%% controller
mdl.tau = @(M) Contact(M,[xs,zs],rs);

%%
mdl.q0(1)   = 0;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation

%[rig] = setupRig(M,L,Modes);
figure

sphere_gmodel.bake.render();hold on;

for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)
    %rig = rig.computeFK(mdl.Log.q(ii,:));
    %rig = rig.update();
    p = shp.FK(mdl.Log.q(ii,:));
    %g = shp.string(mdl.Log.q(ii,:));
    h = plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(45,-15)
    drawnow();
    delete(h)
end
h = plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));

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
function [tau,error] = Contact(mdl,pC,rC)
n = size(mdl.Log.q,1);
t = mdl.Log.t;

error=[];
% Hunt-crossley param
stiffness = 3e-4;
damping_base = 2e-7;

% spikes pos wrt backbone
spike = [eye(3) [0;0;-5]; zeros(1,3) 1];

% Init
tau        = zeros(n,1);

[g,J] = mdl.Shapes.string(mdl.Log.q);

spikes = pagemtimes(g,spike);

positionBackbone = reshape(g(1:3,4,:),3,[]);
positionSpike = reshape(spikes(1:3,4,:),3,[]);

for i = 1:size(positionBackbone,2)
    pA = positionBackbone([1,3],i);
    pB = positionSpike([1,3],i);
    [s,pContact] = ComputeIntersection(pA.',pB.',pC,rC);
    if s && ~isempty(pContact)
        depth = norm(pB-pContact);
        vContactBB = pA - pContact;
        unitVec = -vContactBB/norm(vContactBB); % unit vector pointing from contact to backbone
        unitVec = [unitVec(1);0;unitVec(2)]; % Make it 3D
        body_velo_twist = J(:,:,i)*mdl.Log.dq;
        spatial_velo_twist = Admap(g(1:3,1:3,i),g(1:3,4,i))*body_velo_twist;
        spatial_velo = isomse3(spatial_velo_twist)*[spikes(1:3,4,i);0];
        dd = spatial_velo(1:3).'*unitVec ; % velocity on the spike direction
        spatial_damping_force = -damping_base*(depth^1.1)*dd*unitVec;
        spatial_stiffness_force=(depth)*stiffness*unitVec;
        
        body_force = [zeros(3,1);spikes(1:3,1:3,i).'*(spatial_stiffness_force+spatial_damping_force)];
        tau = tau + J(:,:,i).' * body_force;
    end
end



% for i = 1:size(position,2)
%     if dist(i) <= 0
%         body_velo_twist = J(:,:,i)*mdl.Log.dq;
%         spatial_velo_twist = Admap(g(1:3,1:3,i),g(1:3,4,i))*body_velo_twist;
%         spatial_velo = isomse3(spatial_velo_twist)*[g(1:3,4,i);0];
%         dd = N(i,:)*spatial_velo(1:3);
%         spatial_damping_force = -damping_base*((-dist(i))^1.1)*dd*N(i,:).';
%         spatial_stiffness_force=(-dist(i)^3)*stiffness*N(i,:).';
%         
%         body_force = [zeros(3,1);g(1:3,1:3,i).'*(spatial_stiffness_force+spatial_damping_force)];
%         tau = tau + J(:,:,i).' * body_force;
%     end
% end
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


%% Line intersect circle
function [flag, points] = ComputeIntersection(pA,pB,pC,rC)
    % project C onto AB
    vAB = pB-pA;
    points = [];
    flag = 0;
    Cproj = ProjectionOnLine([pA; pB],pC).';
    d = norm(pC-Cproj);
    if d > rC
        flag = 0; return
    else
        uAB = vAB/norm(vAB);
        tmp = sqrt(rC^2-d^2);
        if tmp == 0
            points = Cproj;
            flag = 1;
        else
            point1 = Cproj + tmp*uAB;
            point2 = Cproj - tmp*uAB;
            if dot(point1-pA,point1-pB) < 0
                points = [points,point1];
                flag = 1;
            end
            if dot(point2-pA,point2-pB) < 0
                points = [points,point2];
                flag = 1;
            end
        end
    end
                
end


function [ProjPoint] = ProjectionOnLine(vector, q)
% vector = [pA; pB]; 
% q = [pC];
% pA, pB, pC row vectors
p0 = vector(1,:);
p1 = vector(2,:);
a = [-q(1)*(p1(1)-p0(1)) - q(2)*(p1(2)-p0(2)); ...
    -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))]; 
b = [p1(1) - p0(1), p1(2) - p0(2);...
    p0(2) - p1(2), p1(1) - p0(1)];
ProjPoint = -(b\a);
end