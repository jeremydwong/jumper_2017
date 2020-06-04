
function [out1,out2,out3]=odeMuscleJumper(t,state,P)
% function [sol,statep,outstruct]=odeMuscleJumper(t,state,P)
% statep,sol,outstruct
% compute state derivatives of humanoid jumper.
% that is the muscle velocities, and
% equations of motion of inverted 4-pendulum.
% % % %
% INPUTS:   time (s)
%           state (vector)
%           P (struct)
% time
% state=[
%     angles                        [1:4]    radians
%     ddt(angles)                   [5:8]    radians/s
%     footxy                        [9:10]   m
%     d/dt(footxy)                   [11:12] m/s
%     d/dt(relative muscle length)   [13:18] dimensionless [0:1]
%     d/dt(gamma);                   [19:24] dimensionless [0:1]
% ]
%
% P: parameter structure,
%     as organized by a=getJumperParams;P=getJumperStruct(a);
% % % %
% OUTPUTS:  (sol)ution vector
%           statep
%           outstruct: useful internal states
% n_m = max(fse-fpe,zeros(size(fmax)));
% mom_m = mustor;
% ese;  tendon excursion            m
% fse;  force in the tendon         N
% rloi; length origin-insertion     m
% q;    active state of muscle      [0:1]
% stim; stim level of muscle.       [0:1]
% % % %
% %%%%% Brief Description %%%%%
% Hill-type computation of derivative of muscle length:
% we consider the muscle-tendon junction. at that point, assuming
% that this point has no mass or that m*a is relatively small compared to
% sum(F), then we assume
% sum(F) = 0.
% F = F_CE + F_SE + F_PE = 0.
% SE and PE are passive springs for which we only need to know lengths.
% so F_CE = -(F_SE+F_PE).
% F_CE is more complex; F_CE = q*FL*FV.
% But we know
% q is a function of stim and length,
% FL is a function of length,
% then we can solve for FV.
% Now, FV is a monotonic function. Thus, we can invert it to get
% muscle-velocity.
% NOTES:
% 0: DK Pai on lumped muscle models(2010).
% 1: the formulation q*FL*FV might be inappropriate (2013 paper Yeo).
% 2: assuming the muscles have no mass results in erroneous equations of
% motion (2010 Pai).

% model STIM as either bang-bang (off-on) or linear-ramp with slope 5.
STIM_BANGBANG = 1;
STIM_SLOPE5 = 2; %linear activation; slope = 5.
typeStim = STIM_BANGBANG;

state=state(:);
% set the muscle stimulation.
cur_stim = [];
if typeStim ==STIM_BANGBANG
    cur_stim=P.istim+(t>P.tstart).*(1-P.istim);
elseif typeStim == STIM_SLOPE5
    addstim = (t>P.tstart).*(t-P.tstart)*5;
    cur_stim=P.istim + addstim;
    cur_stim=cur_stim - (cur_stim>1).*(cur_stim-1);
end
cur_stim = cur_stim(:)';
% /set the muscle stimulation.

% sort out parameters and inputs
nmus=length(P.m.rlsenul);

% Coefficients for ddt(gamma)
% gammap=rm.*(c.*stim-gamma);
rm=P.m.rm;
c= P.m.c;

% Coefficients rho equation for length-dependent calcium sensitivity.
% rho=getal.*((rldak.^ints)-1.0)./(((rldak./lcerel).^ints)-1.0);
rldak= P.m.rldak;
ints= P.m.ints;
getal= P.m.getal;

% Coefficients for computation of active state.
% q=(q0+(rho.*gamma).^exphat)./(1.0+((rho.*gamma).^exphat));
q0= P.m.q0;
exphat= P.m.exphat;

% Coefficients for the force-velocity parameterization
arel= P.m.arel;
brel= P.m.brel;
slopfac= P.m.slopfac;
fasympt= P.m.fasympt;
sloplin= P.m.sloplin;

% Coefficients that determine muscle force per muscle-type.
rlsenul= P.m.rlsenul;       % SE
strainmax= P.m.strainmax;   % SE
rlpenul= P.m.rlpenul;       % PE
stresmax= P.m.stresmax;     % PE
fmax=  P.m.fmax;            % CE csa-based fmax of each muscle.
rlceopt= P.m.rlceopt;       % CE optimum fibre length.
rmapar=  P.m.rmapar;        % CE
width= P.m.width;           % CE
c1= -1.0./(width.^2);       % CE computed
c2= -2.0*c1;                % CE computed
c3=  1.0+c1;                % CE computed

kse=fmax./((strainmax.*rlsenul).^2);
kpe=(stresmax.*fmax)./(((width+1.0-rlpenul).*rlceopt).^2);
kpe = 0;% here we intentionally set PE force to 0.
% % FOUR EXTRAS! G, AIR, CSTIM.
g=P.sim.g;              % gravity
air=P.sim.air;          % in the air
stim=cur_stim(:)';    % current stimulation.
nseg = P.sk.nseg;       % number of segments. remove looping over 4.

% Cofficients for skeletal dynamics
l=P.sk.l;
d=P.sk.d;
p=l-d;
mass= P.sk.mass;
j=P.sk.j;
m=P.sk.mass;
fi=state(1:nseg);
fip=state(nseg+1:2*nseg);
col=2*nseg+4;
lcerel=state(col+1:col+nmus);
gamma=state(col+nmus+1:col+2*nmus);

% Ensure row vectors.
fi = fi(:)';
fip = fip(:)';
lcerel=lcerel(:)';
gamma = gamma(:)';
fijo= [fi(1) diff(fi)];

% Compute the length of the muscle-tendon complex and moment arm
rloi=zeros(1,nmus);
rmomarm=zeros(nseg,nmus);
for iseg=1:nseg
    start=(iseg-1)*3+1;
    rloi=rloi+rmapar(start,:)+rmapar(start+1,:).*fijo(iseg)+rmapar(start+2,:).*fijo(iseg)^2;
    rmomarm(iseg,:)=rmapar(start+1,:)+2*rmapar(start+2,:)*fijo(iseg);
end

% Compute gamma (free calcium), rho (LDCS), q (active state).
gammap=rm.*(c.*stim-gamma);
rho=getal.*((rldak.^ints)-1.0)./(((rldak./lcerel).^ints)-1.0);
q=(q0+(rho.*gamma).^exphat)./(1.0+((rho.*gamma).^exphat));

%
ese=rloi-rlsenul-lcerel.*rlceopt;
fse=(ese>0).*kse.*ese.^2;
epe=(lcerel-rlpenul).*rlceopt;
fpe=(epe>0).*kpe.*epe.^2;

fcerel=max((fse-fpe)./fmax,zeros(size(fmax)));
flenrel=max(c1.*lcerel.^2+c2.*lcerel+c3,ones(size(fmax))*1e-5);

if epe>0
    fprintf('WARNING PE (parallel elastic) forces applied.\n');
end

%calculate vfact and adarel
vfact=min(q*3.333333,ones(size(q)));
adarel=(lcerel<=1).*arel+(lcerel>1).*flenrel.*arel;

%	calculation of ce-velocity
for imus=1:nmus
    if fcerel(imus)/q(imus)<=flenrel(imus) % concentric
        
        %       % this is the hill 1938 equation, scaled for stimulation
        vcerel(imus)=-vfact(imus) * (brel(imus)*(flenrel(imus)+adarel(imus))/(fcerel(imus)/q(imus)+adarel(imus))-brel(imus));
        
    else        % eccentric
        p2=-flenrel(imus)*fasympt(imus);
        p1=(-1/slopfac(imus))*vfact(imus)*brel(imus)*(flenrel(imus)+p2)^2/(flenrel(imus)+adarel(imus));
        p3=-p1/(flenrel(imus)+p2);
        if(fcerel(imus)/q(imus)<=-sqrt(-p1/(vfact(imus)*sloplin(imus)))-p2)
            vcerel(imus)=p1/(fcerel(imus)/q(imus)+p2)+p3;                % eccentric
        else
            vcerel(imus)=vfact(imus)*sloplin(imus)*(fcerel(imus)/q(imus)+p2)+p3+2.0*sqrt(-p1*vfact(imus)*sloplin(imus)); % /*linear eccentric*/
        end
    end
end

%set tor
tor_m=(ones(nseg,1)*fse).*rmomarm;
tor=-sum(tor_m');

dc=d.*cos(fi);ds=d.*sin(fi);
pc=p.*cos(fi);ps=p.*sin(fi);
lc=l.*cos(fi);ls=l.*sin(fi);
m=mass;

if ~isfield(P,'phiground') || ~isfield (P,'gstiff') || ~isfield(P,'gdamp')
    P.phiground = 2.743;
    P.gstiff = -30000;
    P.gdamp = -300;
end;

%% flight phase torques, stiffness and damping
doPassive = 1;

if doPassive && air
    if isfield(P,'kjo')
        kStiffness = P.sk.kjo;
        jPassive = P.sk.djo;
        fiRange  = P.sk.fiRange;
        kPassive = (fijo<fiRange(:,1)).*(fiRange(:,1)-fijo).*kStiffness + ...
            (fijo>fiRange(:,2)).*(fiRange(:,2)-fijo).*kStiffness;
        jDamp = state(5:8).*jPassive;
    else
        kStiffness = [0;70;100;100];
        jPassive = [0;20;40;40];        
        fiRange = [0,pi;
            -pi,0;
            0,pi;
            -pi,0];        
        ddt_phi_rel = [0;%we do not damp the toe.
         state(6)-state(5);
         state(7)-state(6);
         state(8)-state(7);];
            kPassive=(fijo<fiRange(:,1)).*(fiRange(:,1)-fijo).*kStiffness + ...
            (fijo>fiRange(:,2)).*(fiRange(:,2)-fijo).*kStiffness;
        jDamp = -ddt_phi_rel.*jPassive;
        %     j_damp = zeros(4,1);
        %     fiRange=[0,2.7438;
        %         -2.2,-0.1;
        %         .2,2.7;
        %         -3,-.1];
%        else
%         kPassive = zeros(4,1);
%         jDamp = zeros(4,1);
    end
else
    kPassive = zeros(4,1);
    jDamp = zeros(4,1);    
end
    %% /flight
    % The following block of equations are N-E 4-pendulum EOM, [X,Y]=f(fi)
    % A*[F1x,F2x,F3x,F4x,F1y,F2y,F3y,F4y,fidp1,fidp2,fidp3,fidp4,xbasedp,ybasedp]=b
    % where A is square, 14*14
    
    %  F1x    F2x    F3x    F4x    F1y    F2y    F3y     F4y      fidp1     fidp2     fidp3       fidp4    xbasedp  ybasedp
    %C  0      1      2      3      4      5      6      7      8           9           10         11        12      13
    %ML 1      2      3      4      5      6      7      8      9          10          11          12        13      14
    A=[ 1     -1      0      0      0      0      0      0    m(1)*ds(1)    0           0           0       -m(1)     0   ;
        0      1     -1      0      0      0      0      0    m(2)*ls(1)  m(2)*ds(2)    0           0       -m(2)     0   ;
        0      0      1     -1      0      0      0      0    m(3)*ls(1)  m(3)*ls(2)  m(3)*ds(3)    0       -m(3)     0   ;
        0      0      0      1      0      0      0      0    m(4)*ls(1)  m(4)*ls(2)  m(4)*ls(3) m(4)*ds(4) -m(4)     0   ;
        0      0      0      0      1     -1      0      0   -m(1)*dc(1)    0           0           0         0      -m(1);
        0      0      0      0      0      1     -1      0   -m(2)*lc(1) -m(2)*dc(2)    0           0         0      -m(2);
        0      0      0      0      0      0      1     -1   -m(3)*lc(1) -m(3)*lc(2) -m(3)*dc(3)    0         0      -m(3);
        0      0      0      0      0      0      0      1   -m(4)*lc(1) -m(4)*lc(2) -m(4)*lc(3) -m(4)*dc(4)  0      -m(4);
        ds(1)  ps(1)  0      0    -dc(1) -pc(1)   0      0    -j(1)         0           0           0         0       0   ;
        0   ds(2)  ps(2)     0      0    -dc(2) -pc(2)   0      0         -j(2)         0           0         0       0   ;
        0      0   ds(3)  ps(3)     0      0    -dc(3) -pc(3)   0           0         -j(3)         0         0       0   ;
        0      0      0   ds(4)     0      0      0    -dc(4)   0           0           0         -j(4)       0       0   ;
        air    0      0      0      0      0      0      0      0           0           0           0        1-air    0   ;
        0      0      0      0     air     0      0      0      0           0           0           0         0      1-air;
        ];    
    
    
      % 
      %torground = -1* ( fi(1) > P.phiground )*(fi(1)-P.phiground)*P.gstiff - (fi(1) > P.phiground )*fip(1)*P.gdamp;
      torground = 0;
    
    b=[m(1)*(-dc(1)*fip(1)^2);
        m(2)*(-lc(1)*fip(1)^2-dc(2)*fip(2)^2);
        m(3)*(-lc(1)*fip(1)^2-lc(2)*fip(2)^2-dc(3)*fip(3)^2);
        m(4)*(-lc(1)*fip(1)^2-lc(2)*fip(2)^2-lc(3)*fip(3)^2-dc(4)*fip(4)^2);
        m(1)*(-ds(1)*fip(1)^2)-m(1)*g;
        m(2)*(-ls(1)*fip(1)^2-ds(2)*fip(2)^2)-m(2)*g;
        m(3)*(-ls(1)*fip(1)^2-ls(2)*fip(2)^2-ds(3)*fip(3)^2)-m(3)*g;
        m(4)*(-ls(1)*fip(1)^2-ls(2)*fip(2)^2-ls(3)*fip(3)^2-ds(4)*fip(4)^2)-m(4)*g;
        -tor(1)+tor(2)+torground - kPassive(1)+kPassive(2)-jDamp(1)+jDamp(2);
        -tor(2)+tor(3)-kPassive(2)+kPassive(3)-jDamp(2)+jDamp(3);
        -tor(3)+tor(4)-kPassive(3)+kPassive(4)-jDamp(3)+jDamp(4);
        -tor(4)-kPassive(4)-jDamp(4); %+tor(5);
        0;
        0;
        ];
    
    sol=A\b;
    %%%states dphi  ddphi       dbase         ddbase     dlcerel    dgamma
    statep=[fip(:)' sol(9:12)' state(11:12)' sol(13:14)' vcerel(:)' gammap(:)'];
    statep=statep(:); %we output as a column vector.
    outstruct.n_m = max(fse-fpe,zeros(size(fmax)));
    outstruct.tor_m = tor_m;
    outstruct.ese = ese;
    outstruct.fse = fse;
    outstruct.rloi = rloi;
    outstruct.q = q;
    outstruct.stim = stim(:);
    outstruct.glenrel = flenrel;
    outstruct.glenrelq = flenrel.*q;
    outstruct.tor = tor;
    %statep=zeros(size(state));
    
    
    out1 = statep;
    out2 = sol;
    out3 = outstruct;
