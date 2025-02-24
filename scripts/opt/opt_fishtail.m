clr;
%% generate mesh from sdf
sdf = @(x) Trapzoidal(x,15,20,15);

msh = Mesh(sdf,'BdBox',[0,15,0,20],'Quads',[25 50]);
msh = msh.generate().show();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.2,'Penal',4,'FilterRadius',3,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',70,'ChangeMax',0.02);

%% set spatial settings
%fem = fem.set('Repeat',[1 1 1 1],'Periodic',[1, 0]);

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Support',id,[1,0]);
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);

id = fem.FindElements('Location',[7,7],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[7,7],.15);

%% material
fem.Material = TPU90;

%% solving
fem.optimize();
fem.show('ISO');

function Dist = Trapzoidal(P,W,H,E)
R1 = dRectangle(P,0,W,0,H);
L1 = dLine(P,W,0,E,H);
Dist = dIntersect(R1,L1);
end