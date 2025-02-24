function colmap = bluesea(n)
colmap = ...
    [0.1020    0.1373    0.4941
    0.1075    0.1736    0.5041
    0.1127    0.2037    0.5138
    0.1178    0.2296    0.5234
    0.1226    0.2523    0.5327
    0.1273    0.2731    0.5418
    0.1318    0.2921    0.5507
    0.1361    0.3097    0.5594
    0.1404    0.3263    0.5680
    0.1444    0.3418    0.5764
    0.1484    0.3566    0.5846
    0.1523    0.3707    0.5927
    0.1561    0.3841    0.6007
    0.1598    0.3970    0.6085
    0.1634    0.4094    0.6162
    0.1669    0.4213    0.6238
    0.1703    0.4329    0.6313
    0.1737    0.4440    0.6386
    0.1770    0.4548    0.6459
    0.1802    0.4654    0.6530
    0.1834    0.4756    0.6600
    0.1865    0.4856    0.6670
    0.1896    0.4953    0.6738
    0.1926    0.5047    0.6806
    0.1956    0.5140    0.6873
    0.1985    0.5230    0.6939
    0.2013    0.5319    0.7004
    0.2041    0.5406    0.7068
    0.2069    0.5490    0.7132
    0.2097    0.5574    0.7195
    0.2124    0.5656    0.7257
    0.2150    0.5736    0.7318
    0.2176    0.5815    0.7379
    0.2202    0.5892    0.7439
    0.2228    0.5969    0.7499
    0.2253    0.6044    0.7558
    0.2278    0.6117    0.7616
    0.2303    0.6190    0.7674
    0.2327    0.6262    0.7731
    0.2351    0.6332    0.7788
    0.2375    0.6402    0.7844
    0.2398    0.6471    0.7900
    0.2422    0.6538    0.7955
    0.2445    0.6605    0.8009
    0.2467    0.6671    0.8063
    0.2490    0.6736    0.8117
    0.2512    0.6800    0.8170
    0.2534    0.6864    0.8223
    0.2556    0.6927    0.8275
    0.2577    0.6989    0.8327
    0.2741    0.7047    0.8348
    0.3025    0.7103    0.8338
    0.3279    0.7158    0.8328
    0.3512    0.7212    0.8318
    0.3727    0.7266    0.8308
    0.3928    0.7319    0.8298
    0.4118    0.7372    0.8288
    0.4297    0.7424    0.8278
    0.4467    0.7476    0.8268
    0.4630    0.7527    0.8258
    0.4785    0.7578    0.8248
    0.4935    0.7629    0.8238
    0.5079    0.7679    0.8228
    0.5218    0.7728    0.8218
    0.5353    0.7778    0.8208
    0.5483    0.7826    0.8198
    0.5610    0.7875    0.8188
    0.5733    0.7923    0.8177
    0.5852    0.7970    0.8167
    0.5969    0.8018    0.8157
    0.6083    0.8065    0.8147
    0.6194    0.8111    0.8137
    0.6303    0.8157    0.8126
    0.6409    0.8203    0.8116
    0.6513    0.8249    0.8106
    0.6615    0.8294    0.8096
    0.6715    0.8339    0.8085
    0.6813    0.8383    0.8075
    0.6909    0.8427    0.8065
    0.7004    0.8471    0.8054
    0.7097    0.8515    0.8044
    0.7188    0.8558    0.8033
    0.7278    0.8601    0.8023
    0.7367    0.8644    0.8013
    0.7454    0.8687    0.8002
    0.7539    0.8729    0.7992
    0.7624    0.8771    0.7981
    0.7707    0.8813    0.7971
    0.7789    0.8854    0.7960
    0.7870    0.8895    0.7950
    0.7950    0.8936    0.7939
    0.8029    0.8977    0.7928
    0.8107    0.9017    0.7918
    0.8184    0.9058    0.7907
    0.8260    0.9098    0.7897
    0.8335    0.9137    0.7886
    0.8410    0.9177    0.7875
    0.8483    0.9216    0.7865
    0.8556    0.9255    0.7854
    0.8627    0.9294    0.7843];

if nargin>0
    colmap = clrmapping(colmap,n);
end

end

function colmap = clrmapping(colmap,arg)

if arg > 1
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,arg);
    colmap = interp1(x,colmap,xq);
elseif arg < 0
    if abs(arg) == 1, arg = 100; end
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,abs(arg));
    colmap = interp1(x,colmap,xq);
    colmap = flipud(colmap);
elseif arg == 0
    colmap = [colmap;flipud(colmap)];
end

end

