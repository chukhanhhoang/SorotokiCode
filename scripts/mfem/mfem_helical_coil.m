clr;
%% set signed distance function
D = 5.0;
R = 4.0;
H = 15.0;
W = 25;

sdf =  cHelix(D,0,0,R,H,W);

%% generate mesh
mag = Mmesh(sdf,'NElem',1000,'WireThickness',1e-1);
mag = mag.generate();

%% compute inductance
mag = Blender(mag,'Curve',{'PCC+',1e-3,0.0,1.5});
mag.inductance();

%% render
mag.render()
mag.ground();

