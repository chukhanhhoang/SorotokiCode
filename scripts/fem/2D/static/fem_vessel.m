clr;
%% set signed distance function
R = 4;
r = 3;
sdf = @(x) QuarterVessel(x,R,r);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,R,0,R],'NElem',50,'Triangulate',true);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/150);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);

id = fem.FindEdges('EdgeSelect',(R-r)*[cos(pi/6),sin(pi/6)],60);
fem = fem.AddConstraint('Pressure',id,[15*kpa,0]);

% add output nodes
fem = fem.AddConstraint('Output',id{:},[0,0]);

%% assign material
fem.Material = Dragonskin10(100);

%% solving
fem.solve();

function Dist = QuarterVessel(P,R,r)
  C1 = dCircle(P,0,0,R);
  C2 = dCircle(P,0,0,r);
  R1 = dRectangle(P,0,R,0,R);
  
  Dist = dIntersect(dDiff(C1,C2),R1);
end
