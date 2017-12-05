function p = get_jumper_struct(param)
% function p = getJumperStruct(param)
% input: parameter matrix from getJumperParams
% returns parameter structure with two main substructs, m and sk.
% m (23 fields, muscle parameters):
% .nmus = number of muscles
% .names = {'sol','gas','vas','rec','glu','ham'};
% .vals = {'2','23','3','34','4','34'}; moment arms at joints
% 
% The following parameters are associated with the various aspects of the
% hill-type model from the VU. there are conceptually different chunks.
% 1: series-elastic element, SE, the tendon
% (passive quadratic spring)
% 2: a parallel elastic element (passive quadratic spring), PE
% 3:actuating contractile element, CE
% with properties:
% calcium dynamics CD, activation dynamics AD,
% force-length FL, and force-velocity FV.
% 
% CalciumDynamics are a function of Stimulation.
%       d/dt(Calcium)=rm*c*(stim-Calcium)
% Calcium Dynamics parameters (dimensionless, muscle-group independent):
% .rm,     
% .c  
% 
% activeState,'q', Hatze 1981.
% activateState is a function of lcerel(relative muscle length), and [Ca2+]
% rho=getal.*((rldak.^ints)-1.0)./(((rldak./lcerel).^ints)-1.0);
% q=(q0+(rho.*gamma).^exphat)./(1.0+((rho.*gamma).^exphat));
% .q0 minimum active-state of the muscle (i.e. the muscle is always a bit
% on)     
% ActiveState parameters(dimensionless, muscle-group independent):
% .rldak)= 2.9
% .ints= [1] in theory can change getal, rldak, ints to higher order model.
% .getal= 52700
% .exphat=3;
% 
% CE: parameters (some group dependent, as noted):
% .gmax= CE (muscle-group dependent)  
% .rlceopt= relative muscle length
% .width= width of parabola (independent) FL 
% .arel= 0.41; maximum muscle velocity (independent) FV
% .brel= 5.2; relative muscle force (independent) FV
% .slopfac= 2 eccentric slope is slopfac * as steep as concentric FV
% .fasympt= 1.5 * q*gmax is the asymptote line FV
% .sloplin= 200. slope of the line FV
% 
% SE: parameters (independent)
% .rlsenul= SE: resting SE length 
% .straimax= SE: relative strain of the tendon 
% 
% PE: parameters (independent)
% .rlpenul= PE: relative resting length of PE
% .stresmax= PE force is .5*gmax at 1+width. 
% 
% MUSCLE LENGTH PARAMETER MATRIX rmapar:
% .rmapar= matrix of parameters. KvS fits quadratic models of total muscle
% length (loi, length of origin to insertion) as functions of joint angle.
% 
% % % % % % % % % % % % % 
% 
% sk (13 fields, skeleton parameters):
% .jnames = {'toe','ank','kne','hip'};
% .segnames = {'Feet','Shanks','Thighs','HAT'};
% .l= muscle length
% .d= centre of mass from distal joint
% .p= centre of mass from proximal joint (not used)
% .mass=  kg
% .j= moment of inertia
% 
% (to implement flight phase) we implement the following hard-stops:
% .filock= radians; hard-stops of range of motion 
% .kjo= N/m: stiffness of hard-stop forces 
% .djo= N/s/m; damping of hard-stop forces
% 
% Embarrassingly we are mis-using this structure to also house
% initial torques:
% .tor=   param(2,2+1:2+nseg); this need not contain anything
% gravity: 
% .g=param(2,1);
% number of segments:
% .nseg = 4;
 
p=struct;
m = struct;
sk = struct;
mt=param(1,1);
nseg=param(1,2);
segpar=param(4+1:4+mt,1:nseg);

sk.jnames = {'toe','ank','kne','hip'};
sk.segnames = {'Feet','Shanks','Thighs','HAT'};
sk.l=     segpar(1,1:nseg);
sk.d=     segpar(2,1:nseg);
sk.p=sk.l-sk.d;
sk.mass=  segpar(3,1:nseg);
sk.j=     segpar(4,1:nseg);
sk.filock=segpar([5 5+nseg-1],:)';
sk.kjo=   segpar([6 6+nseg-1],:)';
sk.djo=   segpar([7 7+nseg-1],:)';
sk.tor=   param(2,2+1:2+nseg);
sk.g=param(2,1);
sk.nseg = 4;

% % param(1,1) is 12;param(1,3) is 34(index of muspar beginning?)
cm=4+param(1,1);
mt=param(1,3);
nmus=param(1,4);
muspar=param(cm+1:cm+mt,1:nmus);
m.nmus = 6;
m.names = {'sol','gas','vas','rec','glu','ham'};
m.vals = {'2','23','3','34','4','34'};
m.rm=       muspar(1 ,1:nmus);
m.c=        muspar(2 ,1:nmus);
m.q0=       muspar(3 ,1:nmus);
m.rldak=    muspar(4 ,1:nmus);
m.ints=     muspar(5 ,1:nmus);
m.getal=    muspar(6 ,1:nmus);
m.exphat=   muspar(7 ,1:nmus);
m.gmax=     muspar(8 ,1:nmus);
m.rlceopt=  muspar(9 ,1:nmus);
m.width=    muspar(10,1:nmus);
m.arel=     muspar(11,1:nmus);
m.brel=     muspar(12,1:nmus);
m.slopfac=  muspar(13,1:nmus);
m.fasympt=  muspar(14,1:nmus);
m.sloplin=  muspar(15,1:nmus);
m.rlsenul=  muspar(16,1:nmus);
m.straimax= muspar(17,1:nmus);
m.rlpenul=  muspar(18,1:nmus);
m.stresmax= muspar(19,1:nmus);
m.rmapar=   muspar(20:31,1:nmus);

p.m = m;
p.sk = sk;