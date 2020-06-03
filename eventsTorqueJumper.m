function [vnew,isterminal,direction] = eventsTorqueJumper(t,state,P)
% function out = eventsMuscleJumper(t,state,P)
% handle flags. 
if P.i_mode ==1
    [~,sol]=P.ode_handle(t,state,P);
    vnew = sol(5);
    isterminal = 1;
    direction = -1;
else
    vnew = t-P.t_max;
    isterminal = 1;
    direction = 1;
end;
