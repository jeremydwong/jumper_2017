function [statep,sol,outstruct]=odeTorqueJumper(t,state,P)
% function [sol,statep,outstruct]=sol4P(t,state,P)
% compute state derivatives of system.
% system is a torque-driven quadruple-inverted pendulum.
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
% must contain link parameters:
% l, d, mass, j
% and commands:
% tau
% 
% % % %
% OUTPUTS:  (sol)ution vector
%           statep
%           outstruct: useful internal states


l = P.sk.l;           % segment length.
d = P.sk.d;           % the distance from distal to COM. 
p = l-d;              % 
m = P.sk.mass;    %
j =   P.sk.j;        % rotational inertia.

nseg = 4;
state=state(:);
fi = state(1:nseg);fi = fi(:)';
fijo= [fi(1) diff(fi)]; %relative angles.

fip = state(nseg+1:2*nseg);fip = fip(:)';

% [~,i_tor] = min(abs(t-P.t_tor));
% tor = P.tor(i_tor,:);
for iu = 1:nseg
    tor(iu) = ppval(P.U{iu},t);
end
    
if length(tor)==0
    fprintf('error! no torque!\n');
end;

dc = d.*cos(fi);ds=d.*sin(fi);
pc = p.*cos(fi);ps=p.*sin(fi);
lc = l.*cos(fi);ls=l.*sin(fi);

if ~isfield(P,'phiground') || ~isfield (P,'gstiff') || ~isfield(P,'gdamp')
    P.phiground = 2.743;
    P.gstiff = -30000;
    P.gdamp = -300;
end;

%%% flight phase torques, stiffness and damping
doPassive = 1;
torPassive = [0;0;0;0];

air = P.sim.air;
if doPassive && air
    tor = [0;0;0;0];%turn off joint torques once you leave the ground
    fiRange=[0,2.7438;
        -2.2,-0.1;
        .2,2.7;
        -3,-.1];
    kPassive = [0;70;100;100];
    jPassive = [0;20;40;40];    
    fiRange = [0,pi;
        -pi,0;
        0,pi;
        -pi,0];
%     if isfield(P,'kjo')
%         kPassive = P.sk.kjo;
%         jPassive = P.sk.djo;
%         fiRange  = P.sk.fiRange;
%     end;
    
    fijo=fijo(:);
     torPassive= (fijo<fiRange(:,1)).*(fiRange(:,1)-fijo).*kPassive + ...
        (fijo>fiRange(:,2)).*(fiRange(:,2)-fijo).*kPassive + ...
        -fip(:).*jPassive;
% torPassive = 0
    
%     fprintf('t=%.4f',t);
end;
%%% /flight

g = -9.81;
% We aim at the following block of equations:
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



% torground = -1* ( fi(1) > P.phiground )*(fi(1)-P.phiground)*P.gstiff - (fi(1) > P.phiground )*fip(1)*P.gdamp;
torground = -1* ( fi(1) > P.phiground )*(fi(1)-P.phiground)*P.gstiff - (fi(1) > P.phiground )*fip(1)*P.gdamp;

% torground = 0;
b=[m(1)*(-dc(1)*fip(1)^2);
    m(2)*(-lc(1)*fip(1)^2-dc(2)*fip(2)^2);
    m(3)*(-lc(1)*fip(1)^2-lc(2)*fip(2)^2-dc(3)*fip(3)^2);
    m(4)*(-lc(1)*fip(1)^2-lc(2)*fip(2)^2-lc(3)*fip(3)^2-dc(4)*fip(4)^2);
    m(1)*(-ds(1)*fip(1)^2)-m(1)*g;
    m(2)*(-ls(1)*fip(1)^2-ds(2)*fip(2)^2)-m(2)*g;
    m(3)*(-ls(1)*fip(1)^2-ls(2)*fip(2)^2-ds(3)*fip(3)^2)-m(3)*g;
    m(4)*(-ls(1)*fip(1)^2-ls(2)*fip(2)^2-ls(3)*fip(3)^2-ds(4)*fip(4)^2)-m(4)*g;
    -tor(1) + tor(2) + torground - torPassive(1) + torPassive(2);
    -tor(2) + tor(3) - torPassive(2) + torPassive(3);
    -tor(3) + tor(4) - torPassive(3) + torPassive(4);
    -tor(4) - torPassive(4); 
    0;
    0;
    ];

sol=A\b;
if sum(isnan(sol))
    fprintf('t = %.5f',t);
    
    pause;
    
end
%%%states dphi  ddphi       dbase         ddbase     dlcerel    dgamma
statep=[fip(:)' sol(9:12)' state(11:12)' sol(13:14)'];
statep = statep(:);
outstruct.tor = tor;
