% produce single jump from starting initial condition. 
clear all;clc;
plot_lines = {'linewidth',2};
col=get(groot,'DefaultAxesColorOrder');
%
a=load_jumper_params();
P_orig = get_jumper_struct(a);
P = overwrite_params_2017(P_orig);
fi0=[2.5277    0.8295    2.5385    0.7504];
x0=ones(1,6)*.03+rand(1,6)*0.01;
REDO_OPT_STIM = 1;
if REDO_OPT_STIM
    tstim_nonshift = fminsearch(@(x)run4seg_P_optstart(fi0,x,P,0),x0);
else
 %nominal jump stimulation. 
 tstim_nonshift = [0.0217
    0.0208
    0.0204
    0.0221
    0.0010
    0.0120];
end
[height_nonshift,state_nonshift,o_nonshift]=run4seg_P_optstart(fi0,...
    tstim_nonshift,P,0);
tstim_base = tstim_nonshift - min(tstim_nonshift)+.001;
state0 = [fi0,state_nonshift(1,5:end)];
[h_base,state_base,o_base]=run4seg_P(state0,...
                                     o_nonshift.stim(1,:),...
                                     tstim_base,P);

e_base = energy(state_base,o_base,P);


%% animate
[h_vid,state_vid,o_vid]=run4seg_P([fi0,state_nonshift(1,5:end)],...
                                  o_nonshift.stim(1,:),tstim_base,P,1);
animate_jump(o_vid.x,o_vid.y,o_vid.cmx,o_vid.cmy);