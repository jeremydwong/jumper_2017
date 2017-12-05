function param = load_jumper_params
% function param = getJumperParams
% return parameter matrix. used in MFB and KvS code.
load segpar.dat;
[m,nseg]=size(segpar);
param(1,1:2)=[m nseg];
param(4+1:4+m,1:nseg)=segpar;
l=    segpar(1,1:nseg)';
l(4)= 0.5;
d=    segpar(2,1:nseg)';
mass= segpar(3,1:nseg)';
j=    segpar(4,1:nseg)';

load muspar.dat;
[cm,n]=size(param);
[m,nmus]=size(muspar);
param(1,3:4)=[m nmus];
param(cm+1:cm+m,1:nmus)=muspar;
