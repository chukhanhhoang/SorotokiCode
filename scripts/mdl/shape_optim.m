function qd = shape_optim(mdl, object, curve)

shp = mdl.Shapes;
sp = object;
Sd = curve;

e = 0.1;
h = [];
EE = [];
k  = 0;
q = mdl.q0;
while norm(e) > 1.32e-2 && k < 400
    
    % update iteration
    k = k + 1;
    
    % compute Cosserat configuration
    [g, J] = shp.string(q);
    
    % extract positions
    p = reshape(g(1:3,4,:),3,[]).';
    
    % compute closest points
    V = Sd.Node;
    [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % update IK control law
    E(k) = 0;
    dq   = 0;
    
    A = [];
    b = [];
    
    for ii = round(shp.NNode/5):shp.NNode
        
        pd = [XY(ii,1);0;XY(ii,2)];
        Rd = ...%rotz(ii/shp.NNode*2*pi)*...
            -[T(ii,1),B(ii,1),N(ii,1);
              T(ii,3),B(ii,3),N(ii,3);
              T(ii,2),B(ii,2),N(ii,2)];
        
        gd = SE3((Rd),pd);
        
        [dr, dE] = EnergyController(...
            g(:,:,ii), gd,...
            J(:,:,ii),ii/shp.NNode);                     
        
        dq = dq + dr;
        
        E(k) = E(k) + (dE.'*dE);
    end
    
    [Et,R] = shp.tangentPoint(q,...
        [XY(:,1),0*XY(:,1),XY(:,2)]);
    
%     subplot(1,2,1);
%     rig = rig.computeFK(q);
%     rig = rig.update();
%     
%     % setup figure
%     setupFigure(BdBox);
%     
%     hold on;
%     if isempty(h)
%         h = plot3(p(:,1),p(:,2),p(:,3),'b-','LineW',3);
%     else
%         set(h,'XData',p(:,1));
%         set(h,'YData',p(:,2));
%         set(h,'ZData',p(:,3));
%     end
%     
%     subplot(1,2,2);
%     plot(shp.Sigma,R,'LineW',3);
%     axis([0 1 0 0.25]);
%     drawnow;
    
    % compute update state and compute error
    q = q + mdl.G_u*pinv(mdl.G_u)*real(dq);
    e = norm(abs(dq));
     
end

qd = q;

end

function [dq, E] = EnergyController(g,gd,J,sigma) % sigma: 0 to 1
    %
    k1   = 0.003;
    k2   = 10;
    lam1 = 5;
    
    % conditioner
    W  = 1;
    Kp = 5*sigma*diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(W*Xi)*wedge(W*Xi);

    dq = lam1*J.'*Fu;
    
    E = Kp*wedge(Xi);
    
end

function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end