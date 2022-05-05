function qd = nonl_dist_optim(mdl,desired_state)
%     Theta = mdl.Shapes.Ba * mdl.Shapes.get('ThetaEval');
    Theta_3D = pagemtimes(mdl.Shapes.Ba, mdl.Shapes.get('ThetaEval'));
    Theta_cell = num2cell(Theta_3D,[1 2]);
    Theta = vertcat(Theta_cell{:});
    Theta( ~any(Theta,2), : ) = [];
    func = @(x) cost(mdl,Theta,x,desired_state);
    nonlcon = @(x) nonl(mdl,Theta,x,desired_state);
    options = optimoptions('fmincon',... %'SpecifyObjectiveGradient',true,...
        "ConstraintTolerance",2e-8,...
        'StepTolerance',1e-6);
    qd = fmincon(func,desired_state,[],[],[],[],[],[],nonlcon,options);
end

function [f,g] = cost(mdl,Theta,x,desired_state)
    %objective func
    [M_,C_,K_,R_,G_,p_,Phi_,J_,Jt_,Vg_,Kin_] = mdl.getMatrices(x);
    G_annihilator = null(mdl.G_u.').';
    pde_res = G_annihilator*(-G_+K_*x);
    f = pde_res.'*pde_res;
end

function [c,ceq] = nonl(mdl,Theta,x,desired_state)
    c = (x-desired_state).'*(Theta.'*Theta)*(x-desired_state);
    ceq = [];
end