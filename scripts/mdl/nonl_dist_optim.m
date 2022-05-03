function qd = nonl_dist_optim(mdl,obj,desired_state)
    Theta = mdl.Shapes.Ba * mdl.Shapes.get('ThetaEval');
    H = Theta.'*Theta;
    func = @(x) cost(H,x,desired_state);
    nonlcon = @(x) nonl(mdl,x,obj);
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'StepTolerance',1e-10);
    qd = fmincon(func,desired_state,[],[],[],[],[],[],nonlcon,options);
end

function [f,g] = cost(H,x,desired_state)
    %objective func
    f = (x-desired_state).'*H*(x-desired_state);
    if nargout > 1 % gradient required
        g = H*(x-desired_state);
    end
end

function [c,ceq] = nonl(mdl,x,obj)
    [M_,C_,K_,R_,G_,p_,Phi_,J_,Jt_,Vg_,Kin_] = mdl.getMatrices(x);
    

    ceq=[];
end