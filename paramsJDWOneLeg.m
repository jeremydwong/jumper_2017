function P = paramsJDWDangleLeg()
a=paramsLoadMFBsFile2Array();
P_orig = paramsMFBArray2Struct(a);
P = paramsOverwrite2017(P_orig);
massOriginal=P.sk.mass;
P.sk.mass(4) = massOriginal(4)+1/2*sum(massOriginal(1:3));
P.sk.mass(1:3) = 1/2*massOriginal(1:3);%for bobbert 2006.
P.m.fmax = P.m.fmax/2;
