function [height_mus,state,o,tstim,fi22]=mainOptimizeStandardJumper()% optimize muscle stimulations for jump height.
x0=ones(1,6)*.01;
a=getJumperParams;
P = getJumperStruct(a);
P = overwriteparams2017(P);
% initial state of joints
fi22=[2.5277    0.8295    2.5385    0.7504]; %% '22' refers to Bobbert2006 grid.
%%
tstim = zeros(6,1);
REDO_OPT_STIM = 0;
if REDO_OPT_STIM
    tstim = fminsearch(@(x)run4seg_optstart(x,fi22,P,0),x0);
else
    %     [0.1109
%     0.1082
%     0.0980
%     0.1008
%     0.0953
%     0.0879];

% check for a match?
 tstim = [0.0217
    0.0208
    0.0204
    0.0221
    0.0010
    0.0120];
% slightly better jump:1.433
%     tstim =[0.0127
%     0.1269
%     0.0144
%     0.1503
%     0.0508
%     0.0095];

%stims that closely-replicate van soest. in vanSoest1993. 
% t_stim = [0.1109
%     0.1082
%     0.0980
%     0.1008
%     0.0953
%     0.0879];

end;
%get output.
%% simulate the solution.
[height_mus,state,o]=run4seg_optstart(tstim,fi22,P,0);
