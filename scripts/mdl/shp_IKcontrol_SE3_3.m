clr; beep off;
%% settings
L = 100;
M = 4;    % number of modes
N = 50;   % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
    %Y(:,ii) = pcc(x,ii,M);       % piece-wise constant
end

Y = gsogpoly(Y);

%% desired SE3
Sd = sCircle(25,-30,12);

%% soft sobotics shapes
shp = Shapes(Y,[0,M,0,0,0,0],'L0',L);

q  = 1e-4*sort(rand(shp.NDim,1));
q(1) = -0.15;
q(2) = 0.01;
e  = 0.1;

%% solve IK
figure(103); BdBox = [0,75,-50,50];

EE = [];
k  = 0;

Pts = round(linspace(0,N-1,M+1))+1;

while norm(e) > 1e-6 && k < 5000
    clf;
    
    % update iteration
    k = k + 1;
    
    % compute Cosserat configuration
    [g, J] = shp.string(q);
    
    % extract positions
    p  = reshape(g(1:3,4,:),3,[]).';
    
    % compute closest points
    V = Sd.Node;
    [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % plotting
    subplot(1,2,1); hold on;
    %cplane(X,Y,Z);  hold on;
    plot(p(:,1),p(:,3),'k-','LineW',2);
    plot(p(Pts,1),p(Pts,3),'k.','MarkerS',15);
    plot(V(:,1),V(:,2),'k--','MarkerS',10);
    
    % update IK control law
    E(k) = 0;
    dq   = 0;
    
    A = [];
    b = [];
    
    for ii = 1:shp.NNode
        
        pd = [XY(ii,1);0;XY(ii,2)];
        Rd = [T(ii,1),B(ii,1),N(ii,1);
              T(ii,3),B(ii,3),N(ii,3);
              T(ii,2),B(ii,2),N(ii,2)];
        
        gd = SE3((Rd),pd);
        
        [dr, dE] = EnergyController(g(:,:,ii),...
                                    gd,J(:,:,ii),k,ii/shp.NNode);                     
        
        dq = dq + dr;
        
        E(k) = E(k) + (dE.'*dE);
    end
   
    q = q + dq;
    
    % iteration counter
    fprintf(' iteration = %i \n', k);
    
    % compute update state and compute error
    q = q + dq;
    e = norm(abs(dq));
    
    % setup figure
    setupFigure(BdBox);
    
    subplot(1,2,2);
    plot(E,'LineW',3);
%     axis equal;
%     axis([0 100 0 0.05]);
    title('Energy difference');
    drawnow;
   
end

function [dq, E] = EnergyController(g,gd,J,k,i)
    
    k1 = 0.001;
    k2 = i*0.1;
    
    lam1 = 1e-7;
    
    % conditioner
    W  = smoothstep(k/10+0.1);
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(W*Xi)*isomse3(W*Xi);

    dq = lam1*J.'*Fu;
    
    E = Kp*isomse3(Xi);
end

function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end

function setupFigure(B)
    axis equal;
    axis(B);
end
