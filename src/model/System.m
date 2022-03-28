classdef System
    %% Class System
    properties (Access = public)
        Model;
        Object;
        Control;
        
        F_SR;
        F_ob;
        contact_flag=false;
    end
    
    properties (Access = private)
        Contact;
        ContactWrenches;
        ContactPoints;
        
        
    end
    
    %% public methods
    methods
        function sys = System(mdl,obj,contact,control)
            sys.Model = mdl;
            sys.Object = obj;
            sys.Contact = contact;
            sys.Control = control;
            sys.Model.tau = @(mdl) sys.computeTau(mdl);
        end
        
        function [tau,other] = computeTau(System, mdl)
            [F_SR_, ~, flag] = System.Contact(mdl, System.Object);
            [tau,other] = System.Control(mdl, System.Object,flag);
            tau = tau + F_SR_;
        end
        
    end
    
    %% private methods
    methods (Access = private)
        function System = saveContactInfo(System,F_SR_, F_ob_, flag_)
            System.F_SR = F_SR_;
            System.F_ob = F_ob_;
            System.contact_flag = flag_;
        end
        
        
        
    end
    
end