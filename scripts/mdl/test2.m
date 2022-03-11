clr; 
%% 
L = 100;  % length of robot
N = 50;   % number of discrete points on curve
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

%%
mdl = Model(shp,'Tstep',H,'Tsim',1);
mdl.gVec = [0;0;-9.81e3];

%%
mdl.q0(1)   = 0.0;
mdl = mdl.computeEL(mdl.q0);
x=30;y=50;z=-20;
gd = SE3(eye(3),[x;y;z]);
tic
qd = find_qd_from_gL(mdl,gd)
toc
p = shp.FK(qd);

%% 
%% animation
figure;
hold on;
plot3(p(:,1),p(:,2),p(:,3),'LineW',3,'Color',col(1));
scatter3(x,y,z,100);
axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
view(30,30);
drawnow();

function tau = Controller(mdl,r_o)
    qd = forwardbackward(mdl,r_o);
    n = numel(mdl.Log.q);
    t = mdl.Log.t;
    % 
    %tau        = zeros(n,1);

    dV_dq = mdl.Log.EL.G + mdl.Log.EL.K*mdl.Log.q;
    dVd_dq = mdl.Log.EL.K*(mdl.Log.q-qd);
    tau = dV_dq-dVd_dq-4*mdl.Log.EL.M*mdl.Log.dq;
end

function qd = find_qd_from_gL(mdl,gd)
    H = mdl.Log.EL.K;
    func = @(x) cost(H,x);
    nonlcon = @(x) nonl(mdl,x,gd);
    options = optimoptions('fmincon','Display','iter','SpecifyObjectiveGradient',true,'StepTolerance',1e-7);
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

function qd = forwardbackward(mdl,r_o)
    % Init
    n = numel(mdl.q0);
    qd = zeros(n,1);
%     q = mdl.Log.q;
    q = mdl.q0;
    lambda = zeros(7,mdl.Shapes.NNode);
    
    Theta = pagemtimes(mdl.Shapes.Ba,mdl.Shapes.get('ThetaEval'));
    tmp = vertcat(num2cell(Theta, [1 2]));
    Theta_cat = vertcat(tmp{:});
    Ktt = mdl.Shapes.Ktt;
    xi0   = reshape(mdl.Shapes.get('Xi0Eval'),6,[]);
    xi = reshape(pagemtimes(Theta,q),6,[])+xi0;
    ds = mdl.Shapes.ds;
    maxIter = 100;
    rate = 1e-8;
    
    
    for k = 1:maxIter
        % step 1: forward kinematics
        [g,~] = mdl.Shapes.string(q);
        gQ = SE3toQuatR(g);
        % step 2: backward for lambda
        for j = mdl.Shapes.NNode-1:-1:1
            lambda(:,j) = lambda(:,j+1) + ds*dHhatdgQ(gQ(:,j),xi(:,j),lambda(:,j),xi0(:,j),r_o,Ktt).';
        end
        % step 3: update w
        for j = 1:mdl.Shapes.NNode
            xi(:,j) = xi(:,j) + rate*dHhatdxi(gQ(:,j),xi(:,j),lambda(:,j),xi0(:,j),r_o,Ktt).';
        end
        q = Theta_cat\reshape(xi,[],1);
    end
    qd = q;
end

function ret = dHhatdgQ(in1,in2,in3,in4,in5,in6)
% inputs: gQ,xi,lambda,xi0,r_o,Ktt
ret = [(in3(2,:).*in2(1,:))./2.0+(in3(3,:).*in2(2,:))./2.0+(in3(4,:).*in2(3,:))./2.0+in3(7,:).*(in1(2,:).*in2(5,:).*2.0-in1(3,:).*in2(4,:).*2.0)-in3(6,:).*(in1(2,:).*in2(6,:).*2.0-in1(4,:).*in2(4,:).*2.0)+in3(5,:).*(in1(3,:).*in2(6,:).*2.0-in1(4,:).*in2(5,:).*2.0),in3(1,:).*in2(1,:).*(-1.0./2.0)-(in3(3,:).*in2(3,:))./2.0+(in3(4,:).*in2(2,:))./2.0+in3(5,:).*(in1(3,:).*in2(5,:).*2.0+in1(4,:).*in2(6,:).*2.0)-in3(6,:).*(in1(1,:).*in2(6,:).*2.0+in1(2,:).*in2(5,:).*4.0-in1(3,:).*in2(4,:).*2.0)+in3(7,:).*(in1(1,:).*in2(5,:).*2.0-in1(2,:).*in2(6,:).*4.0+in1(4,:).*in2(4,:).*2.0),in3(1,:).*in2(2,:).*(-1.0./2.0)+(in3(2,:).*in2(3,:))./2.0-(in3(4,:).*in2(1,:))./2.0+in3(6,:).*(in1(2,:).*in2(4,:).*2.0+in1(4,:).*in2(6,:).*2.0)+in3(5,:).*(in1(1,:).*in2(6,:).*2.0+in1(2,:).*in2(5,:).*2.0-in1(3,:).*in2(4,:).*4.0)-in3(7,:).*(in1(1,:).*in2(4,:).*2.0+in1(3,:).*in2(6,:).*4.0-in1(4,:).*in2(5,:).*2.0),in3(1,:).*in2(3,:).*(-1.0./2.0)-(in3(2,:).*in2(2,:))./2.0+(in3(3,:).*in2(1,:))./2.0+in3(7,:).*(in1(2,:).*in2(4,:).*2.0+in1(3,:).*in2(5,:).*2.0)-in3(5,:).*(in1(1,:).*in2(5,:).*2.0-in1(2,:).*in2(6,:).*2.0+in1(4,:).*in2(4,:).*4.0)+in3(6,:).*(in1(1,:).*in2(4,:).*2.0+in1(3,:).*in2(6,:).*2.0-in1(4,:).*in2(5,:).*4.0),in1(5,:).*-2.0+in5(1,:).*2.0,in1(6,:).*-2.0+in5(2,:).*2.0,in1(7,:).*-2.0+in5(3,:).*2.0];
end

function ret = dHhatdxi(in1,in2,in3,in4,in5,in6)
% inputs: gQ,xi,lambda,xi0,r_o,Ktt
ret = [(in1(1,:).*in3(2,:))./2.0-(in1(2,:).*in3(1,:))./2.0-(in1(3,:).*in3(4,:))./2.0+(in1(4,:).*in3(3,:))./2.0+in6(1).*(in4(1,:)-in2(1,:)).*2.0+in6(7).*(in4(2,:)-in2(2,:))+in6(13).*(in4(3,:)-in2(3,:))+in6(19).*(in4(4,:)-in2(4,:))+in6(25).*(in4(5,:)-in2(5,:))+in6(2).*(in4(2,:)-in2(2,:))+in6(31).*(in4(6,:)-in2(6,:))+in6(3).*(in4(3,:)-in2(3,:))+in6(4).*(in4(4,:)-in2(4,:))+in6(5).*(in4(5,:)-in2(5,:))+in6(6).*(in4(6,:)-in2(6,:)),(in1(1,:).*in3(3,:))./2.0-(in1(3,:).*in3(1,:))./2.0+(in1(2,:).*in3(4,:))./2.0-(in1(4,:).*in3(2,:))./2.0+in6(7).*(in4(1,:)-in2(1,:))+in6(2).*(in4(1,:)-in2(1,:))+in6(8).*(in4(2,:)-in2(2,:)).*2.0+in6(14).*(in4(3,:)-in2(3,:))+in6(20).*(in4(4,:)-in2(4,:))+in6(26).*(in4(5,:)-in2(5,:))+in6(32).*(in4(6,:)-in2(6,:))+in6(9).*(in4(3,:)-in2(3,:))+in6(10).*(in4(4,:)-in2(4,:))+in6(11).*(in4(5,:)-in2(5,:))+in6(12).*(in4(6,:)-in2(6,:)),(in1(1,:).*in3(4,:))./2.0-(in1(2,:).*in3(3,:))./2.0+(in1(3,:).*in3(2,:))./2.0-(in1(4,:).*in3(1,:))./2.0+in6(13).*(in4(1,:)-in2(1,:))+in6(14).*(in4(2,:)-in2(2,:))+in6(3).*(in4(1,:)-in2(1,:))+in6(9).*(in4(2,:)-in2(2,:))+in6(15).*(in4(3,:)-in2(3,:)).*2.0+in6(21).*(in4(4,:)-in2(4,:))+in6(27).*(in4(5,:)-in2(5,:))+in6(33).*(in4(6,:)-in2(6,:))+in6(16).*(in4(4,:)-in2(4,:))+in6(17).*(in4(5,:)-in2(5,:))+in6(18).*(in4(6,:)-in2(6,:)),in6(19).*(in4(1,:)-in2(1,:))+in6(20).*(in4(2,:)-in2(2,:))+in6(21).*(in4(3,:)-in2(3,:))+in6(4).*(in4(1,:)-in2(1,:))+in6(10).*(in4(2,:)-in2(2,:))+in6(16).*(in4(3,:)-in2(3,:))+in6(22).*(in4(4,:)-in2(4,:)).*2.0+in6(28).*(in4(5,:)-in2(5,:))+in6(34).*(in4(6,:)-in2(6,:))+in6(23).*(in4(5,:)-in2(5,:))+in6(24).*(in4(6,:)-in2(6,:))-in3(5,:).*(in1(3,:).^2.*2.0+in1(4,:).^2.*2.0-1.0)+in3(6,:).*(in1(1,:).*in1(4,:).*2.0+in1(2,:).*in1(3,:).*2.0)-in3(7,:).*(in1(1,:).*in1(3,:).*2.0-in1(2,:).*in1(4,:).*2.0),in6(25).*(in4(1,:)-in2(1,:))+in6(26).*(in4(2,:)-in2(2,:))+in6(27).*(in4(3,:)-in2(3,:))+in6(28).*(in4(4,:)-in2(4,:))+in6(5).*(in4(1,:)-in2(1,:))+in6(11).*(in4(2,:)-in2(2,:))+in6(17).*(in4(3,:)-in2(3,:))+in6(23).*(in4(4,:)-in2(4,:))+in6(29).*(in4(5,:)-in2(5,:)).*2.0+in6(35).*(in4(6,:)-in2(6,:))+in6(30).*(in4(6,:)-in2(6,:))-in3(6,:).*(in1(2,:).^2.*2.0+in1(4,:).^2.*2.0-1.0)-in3(5,:).*(in1(1,:).*in1(4,:).*2.0-in1(2,:).*in1(3,:).*2.0)+in3(7,:).*(in1(1,:).*in1(2,:).*2.0+in1(3,:).*in1(4,:).*2.0),in6(31).*(in4(1,:)-in2(1,:))+in6(32).*(in4(2,:)-in2(2,:))+in6(33).*(in4(3,:)-in2(3,:))+in6(34).*(in4(4,:)-in2(4,:))+in6(6).*(in4(1,:)-in2(1,:))+in6(35).*(in4(5,:)-in2(5,:))+in6(12).*(in4(2,:)-in2(2,:))+in6(18).*(in4(3,:)-in2(3,:))+in6(24).*(in4(4,:)-in2(4,:))+in6(30).*(in4(5,:)-in2(5,:))+in6(36).*(in4(6,:)-in2(6,:)).*2.0-in3(7,:).*(in1(2,:).^2.*2.0+in1(3,:).^2.*2.0-1.0)+in3(5,:).*(in1(1,:).*in1(3,:).*2.0+in1(2,:).*in1(4,:).*2.0)-in3(6,:).*(in1(1,:).*in1(2,:).*2.0-in1(3,:).*in1(4,:).*2.0)];
end




function ret = SE3toQuatR(g)
    if ndims(g) == 2
        ret = [rot2quat(g(1:3,1:3)).';g(1:3,4)];
    elseif ndims(g) ==3
        ret = zeros(7,size(g,3));
        for i = 1:size(g,3)
            ret(:,i) = SE3toQuatR(g(:,:,i));
        end
    end
end

function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end