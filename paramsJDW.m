function P = paramsJDW()
a=paramsLoadMFBsFile2Array();
P_orig = paramsMFBArray2Struct(a);
P = paramsOverwrite2017(P_orig);
P.sk.mass=P.sk.mass;%for bobbert 2006.
