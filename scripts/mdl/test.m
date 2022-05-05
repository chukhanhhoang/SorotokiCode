clr; 
%% 
L = 0.22;  % length of robot
N = 80;   % number of discrete points on curve
M = 6;    % number of modes
H = 1/100; % timesteps
FPS = 10; % animation speed

Modes = [0,M,0,0,0,0];  % pure-XY curvature
%%Object
Sd = sCircle(0.06,-0.06,0.04);  % desired enveloping SDF
obj = sSphere(0.06,0,-0.06,0.04); % offset due to occupance of soft arm
obj_gmodel = Gmodel(obj);
% generate nodal space
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

%%
shp = Shapes(Y,Modes,'L0',L);
shp.E    = 5.00;     % Young's modulus in Mpa
shp.Nu   = 0.49;     % Poisson ratio
shp.Rho  = 100; % Density in kg/m^3
shp.Zeta = 0.1;      % Damping coefficient
shp.Gvec = [0;0;-9.81];
shp = shp.rebuild();
Theta_ = shp.get('ThetaEval');
Theta = pagemtimes(shp.Ba,Theta_);

%%
mdl = cModel(shp,'Tstep',H,'Tsim',6);
% mdl.gVec = [0;0;0-9.81e3];%-9.81e3
mdl.q0(1) = 0;
mdl.q0(3:end) = 2*rand(shp.NDim-2,1);
mdl = mdl.computeEL(mdl.q0);
mdl.G_u = eye(numel(mdl.q0));
mdl.G_u(:,end-2:end) = [];
% find final config
tic
qd = shape_optim(mdl,obj,Sd);
qd = nonl_dist_optim(mdl,qd);
toc
p = shp.FK(qd);




%% controller
% mdl.tau = @(M) Controller_2(M,obj_param);
mdl.tau = @(M) Controller_qd(M,qd);
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
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.q)
    rig = rig.computeFK(mdl.Log.q(ii,:));
    rig = rig.update();
    hold on;
    plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
    scatter3(p(end,1),p(end,2),p(end,3),100,col(2));
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
%     if ii == 1
%         gif('UnderactuatedControl_qd.gif')
%     else
%         gif
%     end
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
        d_ = obj.eval(p(id:end,:));
        d = d_(:,end).^2;
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
function [tau,F_ob] = Controller_2(mdl,sphere)
    F_ob = zeros(6,1);
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    N = mdl.Shapes.NNode;
    sphere_pos = sphere(1:3).';
    sphere_r = sphere(4);
    % 
    %tau        = zeros(n,1);
    % compensate gravity
    dV_dq = mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
%     dV_dq = mdl.Log.EL.K*mdl.Log.q;
    p = mdl.Log.p;
    Phi = mdl.Log.Phi;
    k = 1e-4;
    body_force = zeros(n,1);
    for i = round(N/3):N
        body_force = body_force + mdl.Log.EL.J(:,:,i).'* [zeros(3,1);Phi(:,:,i).'*k*(sphere_pos-p(:,i))];
    end
    
    tau_c = dV_dq -20*mdl.Log.EL.M*mdl.Log.dq + body_force - 5*exp(-1*t) * mdl.Log.EL.K*mdl.Log.q;
    tau_c(end) = 0;
    tau_o = zeros(n,1);
    stiffness = 1e-4;
    damping = 1e-5;
    rs = sphere_r;
    touch = 0;
    for i = 1:N
        vector_from_sphere = p(:,i)-sphere_pos;
        d = norm(vector_from_sphere);
        if d <= rs
            touch = 1;
            body_velo = mdl.Log.EL.J(:,:,i)*mdl.Log.dq;
%             spatial_velo_twist = admap(mdl.Log.Phi(1:3,1:3,i),mdl.Log.p(1:3,i))*body_velo;
%             spatial_velo = wedge(spatial_velo_twist(1:3))*[p(:,i);1];
%             dd = spatial_velo(1:3).'*vector_from_sphere;
%             damp_force = -damping*(rs-d)^1.1*dd*vector_from_sphere/d;
            force = (stiffness*(rs-d)*vector_from_sphere);
            body_force = [zeros(3,1);mdl.Log.Phi(1:3,1:3,i).'*force];
            tau_o = tau_o + mdl.Log.EL.J(:,:,i).' * body_force;
            
%             F_ob = F_ob+ [zeros(3,1); force];
        end
    end
    
%     res = tau-mdl.Log.EL.J(:,:,end).'*[eye(3);zeros(3)]*u
    tau = tau_c+tau_o;
%     tau(3:end) = 0;
end

function [tau,error] = Controller(mdl,sphere)
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    N = mdl.Shapes.NNode;
    sphere_pos = sphere(1:3).';
    sphere_r = sphere(4);
    % 
    %tau        = zeros(n,1);
    error = zeros(6,1);
    % compensate gravity
    dV_dq = -mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
%     dV_dq = mdl.Log.EL.K*mdl.Log.q;
    p = mdl.Log.p;
    Phi = mdl.Log.Phi;
    k = 1e-4;
    body_force = zeros(n,1);
    for i = round(N/3):N
        body_force = body_force + mdl.Log.EL.J(:,:,i).'* [zeros(3,1);Phi(:,:,i).'*k*(sphere_pos-p(:,i))];
    end
    
    tau = dV_dq -10*mdl.Log.EL.M*mdl.Log.dq + body_force;
    
    stiffness = -1e-2;
    damping = 1e-5;
    rs = sphere_r;
    for i = 1:N
        vector_from_sphere = p(:,i)-sphere_pos;
        d = norm(vector_from_sphere);
        if d <= rs
            body_velo = mdl.Log.EL.J(:,:,i)*mdl.Log.dq;
            spatial_velo_twist = Admap(mdl.Log.Phi(1:3,1:3,i),mdl.Log.p(1:3,i))*body_velo;
            spatial_velo = isomse3(spatial_velo_twist)*[p(:,i);1];
            dd = spatial_velo(1:3).'*vector_from_sphere;
            damp_force = -damping*(rs-d)^1.1*dd*vector_from_sphere/d;
            body_force = [zeros(3,1);mdl.Log.Phi(1:3,1:3,i).'*(stiffness*(1-rs/d)*vector_from_sphere+damp_force)];
            tau = tau + mdl.Log.EL.J(:,:,i).' * body_force;
        end
    end
end

%% setup controller
function [tau,error] = Controller_qd(mdl,qd)
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    % 
    %tau        = zeros(n,1);
    error = mdl.Log.q-qd;
    dV_dq = -mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
    dVd_dq = mdl.Log.EL.K*(mdl.Log.q-qd);
    tau = mdl.G_u*pinv(mdl.G_u)*(dV_dq-dVd_dq-3*mdl.Log.EL.M*mdl.Log.dq);
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
rig.g0 = SE3(roty(pi/2)*rotz(pi),zeros(3,1));

rig = rig.render();
end
