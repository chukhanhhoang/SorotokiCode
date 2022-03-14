clr; 
%% 
L = 100;  % length of robot
N = 30;   % number of discrete points on curve
M = 3;    % number of modes
H = 1/125; % timesteps
FPS = 30; % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

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
mdl = Model(shp,'Tstep',H,'Tsim',10);
mdl.gVec = [0;0;-9.81e3];
% mdl = mdl.computeEL(mdl.q0);


%%
mdl = mdl.computeEL(mdl.q0);

%% find final config
tic
qd = find_qd_from_obj(mdl,10);
toc
p = shp.FK(qd(1:M*2));

%% controller
mdl.tau = @(Model) Controller(Model,qd(1:M*2));


mdl = mdl.simulate(); 

%%Object
x = qd(M*2+1);y = qd(M*2+2); z = qd(M*2+3);
obj = sSphere(x,y,z,10);
obj_gmodel = Gmodel(obj);

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.q(:,1:M),'LineW',2);
colororder(col);

%% animation
% figure;
% hold on;
[rig] = setupRig(M,L,Modes);
obj_gmodel.bake().render()
% gif('SimpleGraspControl.gif')
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)
    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    hold on;
    plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
    scatter3(p(end,1),p(end,2),p(end,3),100,col(2));
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
%     gif
end


function qd = find_qd_from_obj(mdl,obj_r)
    p0 = mdl.Shapes.FK(mdl.q0);
    H = mdl.Log.EL.K;
    n = numel(mdl.q0);
    func = @(x) cost(H,x,n);
    nonlcon = @(x) nonl(mdl,x,obj_r,n);
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'StepTolerance',1e-10);
    qd = fmincon(func,[mdl.q0;p0(n,:).'],[],[],[],[],[],[],nonlcon,options);
    
    function [f,g] = cost(H,x,n)
        %objective func
        f = x(1:n).'*H*x(1:n);
        if nargout > 1 % gradient required
            g = [H*x(1:n);zeros(3,1)];
        end
    end

    function [c,ceq] = nonl(mdl,x,obj_r,n)
        robot_state = x(1:n);
        object_pos = x(n+1:end);
        p = mdl.Shapes.FK(robot_state); % positions of points
        id= round(size(p,1)/3); % an third of the robot
        d= vecnorm(bsxfun(@minus,p(2*id+1:end,:).',object_pos));
%         d = d_(:,end);
        c = norm(d-obj_r-2)-1e-3;
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

%% setup controller
function [tau,error] = Controller(mdl,qd)
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    % 
    %tau        = zeros(n,1);
    error = mdl.Log.q-qd;
    dV_dq = mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
    dVd_dq = mdl.Log.EL.K*(mdl.Log.q-qd);
    tau = dV_dq-dVd_dq-8*mdl.Log.EL.M*mdl.Log.dq;
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

