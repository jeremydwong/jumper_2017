function P = getJumperParametersJDW()
a=load_jumper_params();
P_orig = get_jumper_struct(a);
P = overwrite_params_2017(P_orig);
P.sk.mass=P.sk.mass;%for bobbert 2006.
