classdef System
    %% Class System
    properties (Access = public)
        Model;
        Object;
        Control;
    end
    
    properties (Access = private)
        Contact;
        ContactWrenches;
        ContactPoints;
        inContact;
    end
    
    %% public methods
    methods
        function sys = System(mdl,obj,contact,control)
            sys.Model = mdl;
            sys.Object = obj;
            sys.Contact = contact;
            sys.Control = control;
        end
        
        function [flag, F_SR, F_ob] = computeContact(System)
            [flag, F_SR, F_ob] = System.Contact(System.Model, System.Object);
        end
        
        function [Wrenches, Points] = computeControl(System)
            [Wrenches, Points] = System.Control(System.Model, System.Object);
        end
    end
    
    %% private methods
    methods (Access = private)
        
    end
    
end