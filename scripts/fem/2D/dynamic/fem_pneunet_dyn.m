clr;
%% pressure settings
P0  = 20*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;

MinH  = 0.8586*4;
MaxH  = 1.5*4;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/500,'BdBox',[0,120,-80,20],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('SE'),[0,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Contact',@(x) SDF(x,22),[0,0]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
   @(x) P0*sigmoid(x.Time));

%% add output nodes
id  = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = Dragonskin10(75);
fem.Material.Zeta = 1.4;

%% solve
fem.simulate(); 

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,320):numel(t)
%    subplot(1,3,[1,2]);
        fem.set('Node',fem.Log.Node{ii});
        fem.show('Field',fem.Log.Stress{ii});
        axis([-30 120 -100 30]);      
        
%     subplot(1,3,3);
%         cla;
%         K = fem.Log.Kin(1:ii);
%         plot(t(1:ii),K,'LineW',2); hold on;
%         %axis([0 1 -8e-3 7e-3]);
%        
%         U0 = fem.Log.Vg(end) + fem.Log.Psi(end) - fem.Log.Vf(end);
%         U = fem.Log.Psi(1:ii) + fem.Log.Vg(1:ii) - fem.Log.Vf(1:ii) - U0;
%         plot(t(1:ii),U,'LineW',2); 
%         axis([0 t(end) -1 3]);     
% 
%         plot(t(1:ii),K+U,'LineW',2); 
    
    
    background(gitpage);
    drawnow();
    
        
    if ii == 1, gif('fem_pneunet_dyn.gif',...
            'frame',gcf,'nodither','DelayTime',1/40);
    else, gif; 
    end
end

function Dist = SDF(x,R)

sdf = sCircle(30,-50 + R,R);
Dist = sdf.eval(x);
end