clr; 
%% 
L = 100;  % length of robot
N = 80;   % number of discrete points on curve
M = 4;    % number of modes
H = 1/125; % timesteps
FPS = 30; % animation speed

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
shp.Zeta = 0.1;      % Damping coefficient

shp = shp.rebuild();
Theta_ = shp.get('ThetaEval');
Theta = pagemtimes(shp.Ba,Theta_);

%%
mdl = Model(shp,'Tstep',H,'Tsim',2);
mdl.gVec = [0;0;0];
% mdl = mdl.computeEL(mdl.q0);


%%
mdl.q0(1)   = 0.0;
mdl = mdl.computeEL(mdl.q0);

%% find final config
x=10;y=0;z=-10;
gd = SE3(roty(pi/4), [x;y;z]);
% b = isomse3(logmapSE3(shp.get('g0')\gd) - L*isomse3(shp.xia0));
% control_point = L;
% control_point_index = round(control_point/L*N);
% 
% Theta_int = intTheta(Theta,shp.ds);
% 
% A = Theta_int(:,:,control_point_index);
% 
% qd = quadprog( mdl.Log.EL.K ,[],[],[],A,b);
tic
qd = find_qd_from_gL(mdl,gd);
toc
p = shp.FK(qd);

%% controller
mdl.tau = @(M) Controller(M,qd);


mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
% figure;
% hold on;
[rig] = setupRig(M,L,Modes);
% gif('gLtoQdControl.gif')
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)
    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    hold on;
    plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
    scatter3(x,y,z,100,col(2));
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
%     gif
end


function qd = find_qd_from_gL(mdl,gd)
    H = mdl.Log.EL.K;
    func = @(x) cost(H,x);
    nonlcon = @(x) nonl(mdl,x,gd);
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'StepTolerance',1e-7);
    qd = fmincon(func,mdl.q0,[],[],[],[],[],[],nonlcon,options);
    
    function [f,g] = cost(H,x)
        %objective func
        f = x.'*H*x;
        if nargout > 1 % gradient required
            g = H*x;
        end
    end
    
    function g = endeff(mdl,x)
        [g_,~] = mdl.Shapes.string(x);
        g = g_(:,:,end);
    end

    function [c,ceq] = nonl(mdl,x,gd)
        g = endeff(mdl,x);
        c = norm(g(1:3,4)-gd(1:3,4))-1e-3;
        ceq=[];
    end
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
function [tau,error] = Controller(mdl,qd)
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    % 
    %tau        = zeros(n,1);
    error = mdl.Log.q-qd;
    dV_dq = mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
    dVd_dq = mdl.Log.EL.K*(mdl.Log.q-qd);
    tau = dV_dq-dVd_dq-4*mdl.Log.EL.M*mdl.Log.dq;
    tau(end-1:end)=0;
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

