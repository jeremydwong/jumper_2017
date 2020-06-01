function P = getJumperParametersJDWDangleLeg()
a=load_jumper_params();
P_orig = get_jumper_struct(a);
P = overwrite_params_2017(P_orig);
massOriginal=P.sk.mass;
P.sk.mass(4) = massOriginal(4)+1/2*sum(massOriginal(1:3));
P.sk.mass(1:3) = 1/2*massOriginal(1:3);%for bobbert 2006.
P.m.fmax = P.m.fmax/2;