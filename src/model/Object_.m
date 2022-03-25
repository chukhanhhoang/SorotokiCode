classdef Object_
    %% Class System
    properties (Access = public)
        SDF; % SDF class for the outer surface
        m; % mass
        I; % moment of inertia
        pos; % position of the CoM (SE3)
        vel; % velocity of the CoM (body twist)
        r; % radius (temp)
    end
    
    properties (Access = private)
        
    end
    
    methods 
        function obj = Object_(SDF,m,I,pos,vel,r)
            obj.SDF = SDF;
            obj.m = m;
            obj.I = I;
            obj.pos = pos;
            obj.vel = vel;
            obj.r = r;
        end
        
        function d = distance_to_obj(Object_,points)
            % to be filled
        end
    end
    
end