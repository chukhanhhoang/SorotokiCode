classdef cModel
    % Model class with contact
    properties (Access = public)
        Shapes;
        %Inertia;
        Tsim;
        Tstep;

        q, dq, t;
        q0, dq0, Phi0, p0;
        Log;
        constrained_points; constraint_type;planar;
        tau, tau_;
        object;object_center;pts_per_link;
        G_u; % input mapping
        k,npe,r0,rs;       % potential energy power (k/((norm(r-r0)-rs)^n)), n must be odd
    end
    
    properties (Access = private)
        N; S;
        dTaudq, dTauddq;
        
        Xi0;Theta;
        
        MexSolver = true;
        
        ResidualNorm;
        MaxIteration;
        Conv;
        Adaptive;
        
        Linewidth;
        Markersize;
        
        Ba;
        ShpFnc;
    end
    
%--------------------------------------------------------------------------
methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------- cModel Class
function obj = cModel(Shapes,varargin) 
    obj.planar = false;
    obj.Shapes = Shapes;
    obj.Tstep  = 1/60;
    obj.Tsim   = 10;
    obj.Ba     = Shapes.Ba;
    obj.ShpFnc = Shapes.Theta;
    
    obj.MaxIteration = 10;
    obj.ResidualNorm = 0.01;

    G0 = Shapes.get('g0');
    obj.Phi0 = G0(1:3,1:3); 
    obj.p0   = G0(5:7).';    
    obj.q0   = zeros(Shapes.NDim,1) + 1e-3*rand(Shapes.NDim,1);
    obj.dq0  = zeros(Shapes.NDim,1);

    obj.tau  = @(mdl) zeros(Shapes.NDim,1);
    
    obj.Linewidth  = 4;
    obj.Markersize = 25;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end

end
%---------------------------------------------------------------------- get     
function varargout = get(cModel,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = cModel.(varargin{ii});
        end
    else
        varargout = cModel.(varargin);
    end
end   
%---------------------------------------------------------------------- set
function cModel = set(cModel,varargin)
    
    for ii = 1:2:length(varargin)
        if strcmp(varargin{ii},'r0')

        else 
            cModel.(varargin{ii}) = varargin{ii+1};
        end
    end
    
end
%---------------------------------------------------- simulates dyn. system
function cModel = simulate(cModel)
    
cModel.t   = 0:cModel.Tstep:cModel.Tsim;
cModel.q   = []; 
cModel.dq  = [];

cModel.Log    = struct;
cModel.Log.q  = cModel.q0;
cModel.Log.dq = 0.1*cModel.q0;
cModel.Log.EL = struct;

% get evals
cModel.Theta = cModel.Shapes.get('ThetaEval');
cModel.Xi0 = cModel.Shapes.get('Xi0Eval');


if isempty(cModel.dTaudq)
    cModel.Log.t = 0;
    [cModel.dTaudq,cModel.dTauddq] = computeControlJacobians(cModel);
end

[T, X, U, Km, Ue, Ug] = simulateSoftRobot(cModel,...
    [cModel.q0(:); cModel.dq0(:)]);

% extracting data
cModel.Log.t   = T(:);
cModel.Log.q   = X(:,1:cModel.Shapes.NDim,:);
cModel.Log.dq  = X(:,2*cModel.Shapes.NDim+1:cModel.Shapes.NDim);  
cModel.Log.tau = U;

cModel.Log.Kin = Km;
cModel.Log.Vg  = Ug;
cModel.Log.Psi = Ue;

end
%---------------------------------------------------- simulates dyn. system
function cModel = computeEL(cModel,Q,varargin)

if isempty(varargin), dQ = Q*0; end    
cModel.Theta = cModel.Shapes.get('ThetaEval');
cModel.Xi0 = cModel.Shapes.get('Xi0Eval');    
% compute Lagrangian entries
[M_,C_,K_,R_,G_,...
p_,Phi_,J_,Vg_,Kin_] = computeLagrangianFast_mex(...
                                        Q,dQ,... 
                                        cModel.Shapes.ds,...   
                                        cModel.p0,... 
                                        cModel.Phi0,...
                                        cModel.Xi0,... 
                                        cModel.Theta,...
                                        cModel.Shapes.Ba,... 
                                        cModel.Shapes.Ktt,...
                                        cModel.Shapes.Mtt,...     
                                        cModel.Shapes.Zeta,...
                                        cModel.Shapes.Gvec,...
                                        cModel.k,cModel.npe,cModel.r0,cModel.rs);      
                                
% overwrite dynamics
%cModel.Log.t  = T;
cModel.Log.q  = Q;
cModel.Log.dq = dQ;
cModel.Log.p   = p_;
cModel.Log.Phi = Phi_;

cModel.Log.EL.M = M_;
cModel.Log.EL.C = C_;
cModel.Log.EL.R = diag(diag(R_));
cModel.Log.EL.K = diag(diag(K_));
cModel.Log.EL.G = G_;
cModel.Log.EL.J = J_;

cModel.Log.Vg  = Vg_;
cModel.Log.Kin = Kin_;

end
%------------------------------------------------- plot curve configuration
function [P] = show(cModel,Q,col,varargin)
    
    % solve forward kinematics on [0,L0]
    %[P, Np] = computeForwardKinematics(cModel,Q(:));
    
    P = cModel.Shapes.string(Q(:));
    
    if nargin < 2
        col = col(1);
    end
    
    if ~isempty(varargin)
       st = '--'; 
    else
        st = '-';
    end
       
    % plot spatial curve
    plot(P(:,1),P(:,3),'-','linewidth',...
         cModel.Linewidth,'Color',col,'LineStyle',st); hold on;
    % plot interconnection links 
    plot(P([1,end],1),P([1,end],3),'.',...
         'markersize',cModel.Markersize,'Color',col);
end
%--------------------------------------------- compute end-effector pos/vel 
function [p, eta] = endeffector(cModel,Q,dQ)
[p, ~, Jacob] = computeForwardKinematics(cModel,Q(:));
p   = p(end,:).';
eta = Jacob*dQ(:);
end

function [Wt,lambda] = computeConstraints(cModel,Minv,H_)
    % constraint stab constant
    lambda_stab = 40;
    %Jacobian
    J_ = cModel.Log.EL.J;
    Jt_= cModel.Log.EL.Jt;
    % dq
    dQ = cModel.Log.dq;
    
    %compute
    sz2 = cModel.Shapes.NDim;
    n_constraint = length(cModel.constrained_points);
    if n_constraint>0
        if ~cModel.planar
            if cModel.constraint_type == "pinned"
                tmp_1 = [zeros(3,3) eye(3)];
                sz1 = 3;
            else
                sz1 = 6;
                tmp_1 = eye(6);
            end
        else
            if cModel.constraint_type == "pinned"
                sz1 = 2;
                tmp_1 = zeros(2,6);
                tmp_1(1,4) = 1;
                tmp_1(2,6) = 1;
            else
                sz1 = 3;
                tmp_1 = zeros(3,6);
                tmp_1(1,2) = 1;
                tmp_1(2,4) = 1;
                tmp_1(3,6) = 1;
            end
        end
        tmp_2 = zeros(sz1*n_constraint,sz2);
        tmp_3 = zeros(sz1*n_constraint,sz2);

        for i = 1:n_constraint
            tmp_2((i-1)*sz1+1:i*sz1,:) =tmp_1* J_(:,:,cModel.constrained_points(i));
            tmp_3((i-1)*sz1+1:i*sz1,:) =tmp_1*Jt_(:,:,cModel.constrained_points(i));
        end
        Wt = tmp_2;
        wbar = tmp_3*dQ;
        wbar_stab = wbar+lambda_stab*(Wt*dQ);
        lambda = (Wt*Minv*Wt.')\(Wt*Minv*(H_ - cModel.tau_) - wbar_stab);
    end
end

function [points,N,D] = checkContact(cModel)
    P = cModel.pts_per_link;
%     D_ = cModel.object.eval(cModel.Log.p(:,P:P:end).');
%     D = D_(:,end);
%     points = P*find(D<=1.5);
    D_ = cModel.object.eval(cModel.Log.p.');
    D = D_(:,end);
    points = find(D<=1.5);
    N = cModel.Log.p - cModel.object_center;
end

function dQ = resetVelo(cModel)
    J_ = cModel.Log.EL.J;
    sz2 = cModel.Shapes.NDim;
    n_constraint = length(cModel.constrained_points);
    
    % get Wt
    if cModel.constraint_type == "pinned"
        tmp_1 = [zeros(3,3) eye(3)];
        sz1 = 3;
    else
        sz1 = 6;
        tmp_1 = eye(6);
    end
    Wt = zeros(sz1*n_constraint,sz2);

    for i = 1:n_constraint
        Wt((i-1)*sz1+1:i*sz1,:) =tmp_1* J_(:,:,cModel.constrained_points(i));
    end
    
    % quadprog
    M = (cModel.Log.EL.M+cModel.Log.EL.M.')/2;
    x = quadprog(M+1000*eye(size(M)),2*cModel.Log.dq.'*M,[],[],Wt,zeros(sz1*n_constraint,1));
    dQ = cModel.Log.dq + x;
end

function dQ = zeroVelo(cModel)
    dQ = zeros(cModel.Shapes.NDim,1);
end

%%%%%%%%%%%%%%%%%%%%%%%% Get Lagrangian matrices for a static configuration
function [M_,C_,K_,R_,G_,p_,Phi_,J_,Jt_,Vg_,Kin_] = getMatrices(cModel,Q)
    [M_,C_,K_,R_,G_,p_,Phi_,J_,Jt_,Vg_,Kin_]=computeLagrangianFast_mex(...
            Q,zeros(size(Q)),... 
            cModel.Shapes.ds,...   
            cModel.p0,... 
            cModel.Phi0,...
            cModel.Xi0,... 
            cModel.Theta,...
            cModel.Shapes.Ba,... 
            cModel.Shapes.Ktt,...
            cModel.Shapes.Mtt,...     
            cModel.Shapes.Zeta,...
            cModel.Shapes.Gvec,...
            cModel.k,cModel.npe,cModel.r0,cModel.rs);      
end

end
%--------------------------------------------------------------------------
methods (Access = private) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------- implicit time-integration solver
function [Ts, X, U, Kin, Ue, Ug] = simulateSoftRobot(cModel,z0)

Ts  = cModel.t(:);
h   = mean(diff(Ts));
nd  = cModel.Shapes.NDim;

z   = z0(:); 
X   = zeros(length(Ts),nd*2);
U   = zeros(length(Ts),nd);
Kin = zeros(length(Ts),1);
Ue  = zeros(length(Ts),1);
Ug  = zeros(length(Ts),1);

X(1,:)  = z0(:).';
U(1,:)  = zeros(1,nd);

tic;
disp('----------------------------------');
fprintf('* Computing SR dynamics... ');
progress('start')

for ii = 1:length(Ts)-1
    
    % assign states to temp. variable
    dW  = 1e3;
    w_  = z;
    itr = 1;
    
    while (dW > cModel.ResidualNorm && itr <= cModel.MaxIteration)
        
        % compute flow field
        if itr == 1
            [f1, cModel, Mi,reset] = flow(cModel,z,Ts(ii));
            dR = -h*f1;
        else
            [f2, cModel, Mi,~] = flow(cModel,w_,Ts(ii)+h);
            dR = w_ - z - 0.5*h*(f1 + f2);
        end
        
        if isempty(reset)
            % compute hessian
            H = buildHessian(cModel,Mi);

            % hessian update iteration
            dw = stateUpdate(cModel,H,dR);

            % update local state solution
            w_  = w_  + dw;
            itr = itr + 1;

            % compute convergence 
            dW = norm(dw(1:cModel.Shapes.NDim));
        else
            w_  = reset;
            dW = 0;
        end
    end
    
    % update states
    z = w_;
       
    % write output data
    X(ii+1,:) = z(:).';
    U(ii+1,:) = cModel.tau_(:).';
    Kin(ii+1) = 0.5*z(nd+1:2*nd).'*cModel.Log.EL.M*z(nd+1:2*nd);
    Ue(ii+1)  = 0.5*z(1:nd).'*cModel.Log.EL.K*z(1:nd);
    Ug(ii+1)  = cModel.Log.Vg;

    % update progress bar
    inc = round(100*(ii/(length(Ts)))-1);
    progress(inc,100);
end

progress(100,100);
progress('end');
disp('----------------------------------');

tt = toc;
disp(['* Number of elem.  = ', num2str(cModel.Shapes.NNode,3)]);
disp(['* Computation time = ', num2str(tt,3), ' s']);
disp(['* Computation freq = ', num2str(round(1/h)), ' Hz']);
disp(['* Real-time ratio  = ', num2str((Ts(end))/(tt),4), '']);

disp('----------------------------------');

    function [f, cModel, Minv,reset] = flow(cModel,Z,T)
        reset = [];
        n  = cModel.Shapes.NDim;
        Q  = Z(1:n);
        dQ = Z(n+1:2*n);
       
        % compute Lagrangian entries
        cModel = cModel.computeMatrices(Q,dQ,T);
        
        % evaluate control action
        cModel.tau_ = cModel.tau(cModel);
        
        % pre-compute Minverse
        Minv = cModel.Log.EL.M\eye(numel(Q));
        
%         function ret = find_consecutive(A)
%             if size(A,1)~=1
%                 A = A.';
%             end
%             B = [find(diff(A)>1) length(A)];
%             ret = mat2cell(A,1,[B(1) diff(B)]);
%         end
%         
%         % check contact
%         [pts,N,D] = cModel.checkContact();
%         tmp = find_consecutive(pts);
%         ends = [];
%         for i = 1:length(tmp)
%             min_ = min(tmp{i});
%             max_ = max(tmp{i});
% %             mid_ = round((min_+max_)/2);
%             tmp2 = unique([min_,max_]);
%             ends = [ends tmp2];
%         end
% %         two_ends = pts;
%         if length(ends) ~= length(cModel.constrained_points)
% %             cModel.constrained_points = pts;
%             cModel.constrained_points = ends;
%             dQ = cModel.resetVelo();
% %             dQ = cModel.zeroVelo();
%             cModel = cModel.computeMatrices(Q,dQ,T);
%             reset = [Q;dQ];
%         else
%             if any(ends ~= cModel.constrained_points)
%     %             cModel.constrained_points = pts;
%                 cModel.constrained_points = ends;
%                 dQ = cModel.resetVelo();
%     %             dQ = cModel.zeroVelo();
%                 cModel = cModel.computeMatrices(Q,dQ,T);
%                 reset = [Q;dQ];
%             end
%         end
        
        % pre-compute H
        H_ = cModel.Log.EL.C*dQ + cModel.Log.EL.K*Q + cModel.Log.EL.R*dQ - cModel.Log.EL.G;
        % compute constraints
        if ~isempty(cModel.constrained_points)
            [Wt,lambda] = cModel.computeConstraints(Minv,H_);
        else
            Wt = 0;
            lambda = 0;
        end
        
        % flow field
        f = [dQ; ...
             Minv*(cModel.tau_ - H_+ Wt.'*lambda)];
    end
    
end

%------------------------------------- compute Lagrangian matrices and save
function cModel = computeMatrices(cModel,Q,dQ,T)
    if ~cModel.MexSolver
        [M_,C_,K_,R_,G_,...
            p_,Phi_,J_,Vg_,Kin_] = computeLagrangian(cModel,Q,dQ);
    else
        [M_,C_,K_,R_,G_,...
            p_,Phi_,J_,Jt_,Vg_,Kin_] = computeLagrangianFast_mex(...
            Q,dQ,... 
            cModel.Shapes.ds,...   
            cModel.p0,... 
            cModel.Phi0,...
            cModel.Xi0,... 
            cModel.Theta,...
            cModel.Shapes.Ba,... 
            cModel.Shapes.Ktt,...
            cModel.Shapes.Mtt,...     
            cModel.Shapes.Zeta,...
            cModel.Shapes.Gvec,...
            cModel.k,cModel.npe,cModel.r0,cModel.rs);      
    end

    % overwrite dynamics
    cModel.Log.t  = T;
    cModel.Log.q  = Q;
    cModel.Log.dq = dQ;
    cModel.Log.p   = p_;
    cModel.Log.Phi = Phi_;

    cModel.Log.EL.M = M_;
    cModel.Log.EL.C = C_;
    cModel.Log.EL.R = R_;
    cModel.Log.EL.K = K_;
    cModel.Log.EL.G = G_;
    cModel.Log.EL.J = J_;
    cModel.Log.EL.Jt= Jt_;

    cModel.Log.Vg  = Vg_;
    cModel.Log.Kin = Kin_;
end

%--------------------------------------- forwards integration of kinematics
function [pp, id, J] = computeForwardKinematics(cModel,x)
    
% compute total length
ds  = cModel.Shapes.ds;
s    = 0.0;
p_   = cModel.p0;
Phi_ = cModel.Phi0;
J    = zeros(6,numel(x));

pp = p_';

% numerical trapzoid solver over [0,L0]
for ii = 1:cModel.Shapes.NNode-1
    
   s = cModel.Shapes.Sigma(ii);
    
   % first trapzoidal step 
   [K1p,K1Phi,K1J] = ForwardKinematicODE(cModel,s,...
                        x, Phi_, p_);
   % second trapzoid step
   [K2p,K2Phi,K2J] = ForwardKinematicODE(cModel, s + (2/3)*ds,...
                        x, Phi_ + (2/3)*ds*K1Phi, p_ + (2/3)*ds*K1p); 
                    
   p_   = p_   + 0.25*ds*(K1p + 3*K2p);
   Phi_ = Phi_ + 0.25*ds*(K1Phi + 3*K2Phi);
   J    = J    + 0.25*ds*(K1J + 3*K2J);
   
   pp = [pp;p_.'];   
   
end

% transform Jacobian to body frame
J  = adjointSE3inv(Phi_,p_)*J;
id = round(linspace(1,cModel.Shapes.NNode,2));

end
%---------------------------------------------- compute Lagrangian entities 
function [M,C,K,R,G,p,Phi,J,Vg,Kin] = computeLagrangian(cModel,x,dx)

% compute total length
n    = numel(x);
ds   = cModel.Shapes.ds;
p    = cModel.p0;
Phi  = cModel.Phi0;
s    = 0;

% pre-computed Theta and Xi0;
Th = cModel.Theta;
Xi = cModel.Xi0;

Z1 = zeros(6,6+2*(n-1));
Z2 = zeros(n,3*n+1);
Z1(1:3,1:3) = Phi;
Z1(1:3,4)   = p;

if isa(cModel.Shapes.Ktt,'function_handle')
   NLStiff = true; 
else
   NLStiff = false; 
end

for ii = 1:cModel.Shapes.NNode
    
    % first EL-diff eval
    [K1Z1,K1Z2] = LagrangianODEX(cModel,Th(:,:,2*ii-1),Xi(:,1,2*ii-1),...
        x, dx, Z1,NLStiff);
    
    % second EL-diff eval
    [K2Z1,K2Z2] = LagrangianODEX(cModel,Th(:,:,2*ii), Xi(:,1,2*ii),...
        x, dx, Z1 + (2/3)*ds*K1Z1, NLStiff);
    
    % update integrands
    s  = s  + ds;
    Z1 = Z1 + 0.25*ds*(K1Z1 + 3*K2Z1);
    Z2 = Z2 + 0.25*ds*(K1Z2 + 3*K2Z2);

end

% recover the kinematics entities
p   = Z1(1:3,4);
Phi = Z1(1:3,1:3);
B1  = Z1(1:6,5:5+n-1);
J   = Admapinv(Phi,p)*B1;

% recover the dynamics entities
M  = Z2(1:n,1:n);
C  = Z2(1:n,n+1:2*n);
K  = Z2(1:n,2*n+1:3*n);
G  = Z2(1:n,3*n+1);

Vg  = Z1(5,4);
Kin = Z1(6,4);

R = cModel.Shapes.Zeta*K;

end
%-------------------------------------------------- forwards kinematics ODE
function [dp,dPhi,dJ] = ForwardKinematicODE(cModel,s,x,Phi_,p_)

% construct geometric vectors
Theta = cModel.Shapes.Theta(s);
XI    = cModel.Shapes.Ba*Theta*x + cModel.Shapes.xia0;

U     = XI(4:6);
Gamma = XI(1:3);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*skew(Gamma);
A    = adjointSE3(Phi_,p_);
dJ   = A*cModel.Shapes.Ba*Theta;
end

%---------------------------- Lagrangian Matrix-Differential Equation (MDE)
function [dZ1,dZ2] = LagrangianODEX(cModel,Theta_,Xi0_,x,dx,Z1,NLStiff)

n     = numel(x);
p_    = Z1(1:3,4);
Phi_  = Z1(1:3,1:3);
J_    = Z1(1:6,5:5+n-1);
Jt_   = Z1(1:6,6+n-1:6+2*(n-1));

%Theta_ = ThetaEval;%cModel.ShpFnc(s);
XI = cModel.Ba*Theta_*x + Xi0_;

Gamma = XI(1:3);
U     = XI(4:6);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*isomSO3(Gamma);

A   = Admap(Phi_,p_);
Ai  = Admapinv(Phi_,p_);

% build jacobian
Jg  = Ai*J_;
Jgt = Ai*Jt_;

V    = Jg*dx;
adXi = admap(XI);
adV  = admap(V);

BTh = cModel.Shapes.Ba*Theta_;

dJ  = A*BTh;
dJt = A*adV*BTh;
Mtt = cModel.Shapes.Mtt;

if ~NLStiff
    Ktt = cModel.Shapes.Ktt;
else
    Ktt = cModel.Shapes.Ktt(XI);
end

%deta = -adXi*V + BTh*dx;

% compute inertia, coriolis, gravity
dM = (Jg).'*Mtt*Jg;
dC = (Jg).'*((Mtt*adV - adV.'*Mtt)*Jg  + Mtt*Jgt);
dG = (Jg).'*(Ai*Mtt*[0;0;0;0;0;-9.81e3]);

% compute (nonlinear stiffness)
%Ktt = nonlinearStiffnessMat(cModel,s,x);
dK = (BTh).'*Ktt*(BTh);

% compute grav. potential energy
dKe = 0.5*V.'*Mtt*V;
dVg = Mtt(4,4)*p_.'*[0;0;9.81e3];

dZ1                      = zeros(6,6+2*(n-1));
dZ1(1:3,1:3)             = dPhi;
dZ1(1:3,4)               = dp;

dZ1(1:6,5:5+n-1)         = dJ;
dZ1(1:6,6+n-1:6+2*(n-1)) = dJt;

dZ1(5,4)                 = dVg;
dZ1(6,4)                 = dKe;

dZ2 = zeros(n,3*n+1);
dZ2(1:n,1:n)       = dM;
dZ2(1:n,n+1:2*n)   = dC;
dZ2(1:n,2*n+1:3*n) = dK;
dZ2(1:n,3*n+1)     = dG;

end
%--------------------------------------- updated Hessian with dyn. residual
function dr = stateUpdate(cModel, H, dR)
dr = -(-(1/2)*cModel.Tstep*H + eye(size(H,1),size(H,1)))\dR;
end
%---------------------------------------------- compute Hessian approximate
function DF = buildHessian(cModel,varargin)

Minv = varargin{1};
n    = size(cModel.Log.q,1);

DF                      = zeros(2*n,2*n);
DF(1:n,n+1:2*n)         = eye(n);
DF(n+1:2*n,1:n)         = -Minv*(cModel.Log.EL.K - cModel.dTaudq);
DF(n+1:2*n,n+1:2*n)     = -Minv*(cModel.Log.EL.R + cModel.Log.EL.C - cModel.dTauddq);

end
%--------------------------------------------------------- show solver info
function showInformation(cModel)   

fprintf('--------------------------------------------------------------\n');  
fprintf('* Element = %i \n',NodeNum);
fprintf('* Max iteration = %i \n', Fem.MaxIteration);
fprintf('* Solver time horizon = %i \n', Nodel.TimeEnd);
fprintf('* Solver time step    = %i1.1e \n', cModel.TimeStep);
showMaterialInfo(Fem)
fprintf('--------------------------------------------------------------\n');

end
%----------------------------------------- compute nonlinear stiffness mat.
function Ktt = nonlinearStiffnessMat(cModel,s,x)

Jt = cModel.Inertia.Jtt;

E0 = (cModel.E);
G0 = (cModel.E)/(2*(1+cModel.Nu));
  
QE = diag([E0,G0,G0]);
QR = diag([G0,E0,E0]);

Theta1 = cModel.Shapes.Phi(s);
Ba    = cModel.Shapes.Ba;
xi    = Ba*Theta1*x;

%diag([G0*J11,E0*J22,E0*J33,E0*A,G0*A,G0*A]);
Ktt = blkdiag(QR*Jt,QE*eye(3));
Kappa = norm(abs(xi(1:3)));
Gamma = norm(abs(xi(4:6)));

Ktt = Ktt*(1.0 + 1.35*(tanh(-16200*Kappa)^2));

end
%----------------------------------- (pre)-computes the controller jacobian
function [Kt, Dt] = computeControlJacobians(cModel)
n   = length(cModel.q0);
Q0  = cModel.q0(:);    
dQ0 = cModel.dq0(:);    
    
% compute Lagrangian entries
[M_,C_,K_,R_,G_,p_,Phi_,J_] = computeLagrangianFast_mex(...
            Q0,dQ0,... 
            cModel.Shapes.ds,...   
            cModel.p0,... 
            cModel.Phi0,...
            cModel.Xi0,... 
            cModel.Theta,...
            cModel.Shapes.Ba,... 
            cModel.Shapes.Ktt,...
            cModel.Shapes.Mtt,...     
            cModel.Shapes.Zeta,...
            cModel.Shapes.Gvec,...
            cModel.k,cModel.npe,cModel.r0,cModel.rs);
% overwrite dynamics
cModel.q  = Q0; 
cModel.dq = dQ0;
cModel.Log.EL.M  = M_; 
cModel.Log.EL.C  = C_; 
cModel.Log.EL.R  = R_; 
cModel.Log.EL.G  = G_; 
cModel.Log.EL.K  = K_; 
cModel.Log.EL.J  = J_;

cModel.Log.p   = p_; 
cModel.Log.Phi = Phi_; 

cModel.t  = 0;

epsilon = 1e-3;
delta   = epsilon*eye(n);

Tau0 = cModel.tau(cModel);

Kt = zeros(n);
Dt = zeros(n);

% finite difference for tau(q(t),.)
for ii = 1:n
    cModel.q = Q0 + delta(:,ii);
    Tau_ = cModel.tau(cModel);
    Kt(:,ii) = (Tau_ - Tau0)/epsilon;
end

cModel.q  = Q0; 

% finite difference for tau(q(t),.)
for ii = 1:n
    cModel.dq = dQ0 + delta(:,ii);
    Tau_ = cModel.tau(cModel);
    Dt(:,ii) = (Tau_ - Tau0)/epsilon;
end

end

end
end
