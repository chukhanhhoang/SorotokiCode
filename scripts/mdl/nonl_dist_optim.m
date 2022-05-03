function qd = nonl_dist_optim(mdl,desired_state)
%     Theta = mdl.Shapes.Ba * mdl.Shapes.get('ThetaEval');
    
    func = @(x) cost(x,desired_state);
    nonlcon = @(x) nonl(mdl,x);
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,...
        "Algorithm","interior-point",...
        "EnableFeasibilityMode",true,...
        "ConstraintTolerance",1e-3,...
        'StepTolerance',1e-3);
    qd = fmincon(func,desired_state,[],[],[],[],[],[],nonlcon,options);
end

function [f,g] = cost(x,desired_state)
    %objective func
    f = (x-desired_state).'*(x-desired_state);
    if nargout > 1 % gradient required
        g = (x-desired_state);
    end
end

function [c,ceq] = nonl(mdl,x)
    [M_,C_,K_,R_,G_,p_,Phi_,J_,Jt_,Vg_,Kin_] = mdl.getMatrices(x);
    G_annihilator = null(mdl.G_u.').';
    c = [];
    ceq=G_annihilator*(G_+K_*x);
end