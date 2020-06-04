function [height,state,tor_sim,P] = getTorqueModel(state,o,P)
% function [height,tor_sim,P] = getTorqueModel(state,o,P)
% get the torque model from the info. 
t_tor = o.t';
P.tor = o.tor;
P.U = {};
for i =1:4
    P.U{i} = spline(t_tor,P.tor(:,i));
end
% [height,state,tor_sim]=jumpMuscleTorque(state(1,1:12),P,0);
[height,state,tor_sim]=jumpTorque(state(1,1:12),P,0);
