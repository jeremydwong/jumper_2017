function [cost,q,Fcost,Lcecost,Qcost]= cost_FnormFmaxlce(Fs,eqo)
% function [out,Fcost,Lcecost,q]= cost_FnormFmaxlce(Fs,Fs_Max,eqo)
% cost function for computing equilibrium starting positions.
% there are three terms: Fcost, Lcecost, and qcost(1 and 2). 
% Fcost is 
% cost = Fcost+Lcecost+qcost1*qgain+qcost2*qgain;
% 'Lcecost' within legal lengths of muscle (lcrel < 1.4)
% 'qcost[12]' keep within legal q (activate state of muscle) 
fmax = eqo.fmax;
kse = eqo.kse;
rlceopt = eqo.rlceopt;
rloi = eqo.rloi;
rlsenul = eqo.rlsenul;
c1 = eqo.c1;
c2 = eqo.c2;
c3 = eqo.c3;
fse = Fs(:)'; 
ese=sqrt((fse)./kse);
lcerel=(rloi-rlsenul-ese)./rlceopt;
F_isom = lcerel.^2.*c1 + lcerel.*c2 + c3;
Fs = Fs(:);
%1. Fcost is the sum of the forces
Fcost = sum( (Fs ./ (fmax(:) .* F_isom(:) ) ).^2);
Lcecost = sum((lcerel>1.39).*(lcerel-1.39).^2)*2;
qcost1=0;
qcost2=0;

% penalize if q is greater than 1 or less than q0. 
qgain = 100;
q = Fs./(F_isom(:).*fmax(:)); 
qcost1 = sum((q>1).*(q-1).^2);
q0=0.005;
qcost2 = sum((q<q0).*(q-q0).^2);
% end;
Qcost = qcost1*qgain + qcost2*qgain;
cost = Fcost + Lcecost + Qcost;
