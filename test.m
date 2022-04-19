%% Object
xsph = 40; zsph = -50; rsph = 10;
sp=sSphere(50,50,0,10);
sc=sCircle(50,50,10);

[pc,Dc] = distance2curve(sc.Node,[0 0])
[ps,Ds] = distance2curve(sp.Node,[0 0 0])

sp.eval([0 0 0;0 1 0])