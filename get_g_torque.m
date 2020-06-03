function out = get_g_torque(fio,param)
% function out = get_g_torque(fio,p)
% return jount torques from gravity. P must contain
% INPUTS:
% fio - starting position.
% fio is the inertial coordinate system outside the body, 
% angles all relative to rightward x axis. 
% typically, the jumper faces right in this 2D world.
% p - parameter structure, containing
%   .mass (1x4)
%   .l - length of each segment
%   .d - distance to centre of mass from distal endpoint(1x4)
%  [so from toe tip to COM_foot; ankle to COM_shank; knee to COM_thigh; hip to
%  COM_HeadArmTrunk
% OUTPUTS:
% torques at each joint resulting from gravity.
fio = fio(:);
m = param.sk.mass;
g = -9.81;
d = param.sk.d;
l = param.sk.l;
[x,]=xyc4(fio,zeros(size(fio)),zeros(size(fio)),[0;0],[0;0],[0;0],l);
cmxi=[d(1)    0   0     0;
      l(1) d(2)   0     0;
      l(1) l(2) d(3)    0;
      l(1) l(2) l(3) d(4)]*cos(fio);
cmx=m/sum(m)*cmxi;

ctor=[-cmx*sum(m)*g;
    -(m(2:4)/sum(m(2:4)) * (cmxi(2:4)-x(2))) * sum(m(2:4))*g;
    -(m(3:4)/sum(m(3:4)) * (cmxi(3:4)-x(3))) * sum(m(3:4))*g;
    -d(4)                * cos(fio(4))       * m(4)*g];
out = ctor(:);
