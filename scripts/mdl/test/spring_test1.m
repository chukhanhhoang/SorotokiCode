clr; 
%% 
L = 100;  % length of robot
N = 30;   % number of discrete points on curve
M = 3;    % number of modes
H = 1/125; % timesteps
FPS = 30; % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature
%%Object
obj_param = [25,0,-30,12];
obj = sSphere(obj_param(1),obj_param(2),obj_param(3),obj_param(4));
obj_gmodel = Gmodel(obj);
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
mdl = Model(shp,'Tstep',H,'Tsim',20);
mdl.gVec = [0;0;0];%-9.81e3
mdl.q0(1) = 0.125;
% mdl = mdl.computeEL(mdl.q0);


%%
mdl = mdl.computeEL(mdl.q0);

%% find final config
% tic
% qd = find_qd_from_obj(mdl,obj);
% toc
% p = shp.FK(qd);

%% controller
mdl.tau = @(M) Controller(M,obj_param);


mdl = mdl.simulate(); 

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
%     plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
%     scatter3(p(end,1),p(end,2),p(end,3),100,col(2));
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
%     gif
end


function qd = find_qd_from_obj(mdl,obj)
    H = mdl.Log.EL.K;
    func = @(x) cost(H,x);
    nonlcon = @(x) nonl(mdl,x,obj);
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'StepTolerance',1e-10);
    qd = fmincon(func,mdl.q0,[],[],[],[],[],[],nonlcon,options);
    
    function [f,g] = cost(H,x)
        %objective func
        f = x.'*H*x;
        if nargout > 1 % gradient required
            g = H*x;
        end
    end

    function [c,ceq] = nonl(mdl,x,obj)
        p = mdl.Shapes.FK(x); % positions of points
        id= round(size(p,1)/3); % an third of the robot
        d_ = obj.eval(p(2*id:end,:));
        d = d_(:,end);
        c = norm(d-2)-5e-3;
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
function [tau,error] = Controller(mdl,sphere)
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    N = mdl.Shapes.NNode;
    sphere_pos = sphere(1:3).';
    sphere_r = sphere(4);
    % 
    %tau        = zeros(n,1);
    error = [];
    % compensate gravity
    dV_dq = mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
%     dV_dq = mdl.Log.EL.K*mdl.Log.q;
    p = mdl.Log.p;
    Phi = mdl.Log.Phi;
    k = 6e-5;
    body_force = zeros(n,1);
    for i = 2*N/3:N
        body_force = body_force + mdl.Log.EL.J(:,:,i).'* [zeros(3,1);Phi(:,:,i).'*k*(sphere_pos-p(:,i))];
    end
    
    tau = dV_dq -8*mdl.Log.EL.M*mdl.Log.dq + body_force;
    
    stiffness = -1e-2;
    damping = 1e-5;
    rs = sphere_r;
    
    for i = 1:N
        vector_from_sphere = p(:,i)-sphere_pos;
        if norm(vector_from_sphere) <= rs
            body_velo = mdl.Log.EL.J(:,:,i)*mdl.Log.dq;
            spatial_velo = Admap(mdl.Log.Phi(1:3,1:3,i),mdl.Log.p(1:3,i))*body_velo;
            body_force = [zeros(3,1);mdl.Log.Phi(1:3,1:3,i).'*(stiffness*(1-rs/norm(vector_from_sphere))*vector_from_sphere-damping*spatial_velo(4:end))];
            tau = tau + mdl.Log.EL.J(:,:,i).' * body_force;
        end
    end
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

