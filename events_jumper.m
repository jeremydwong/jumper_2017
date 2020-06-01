function [vnew,isterminal,direction] = events_jumper(t,state,P)
% function out = events_jumper(t,state,P)
% handle flags. 
if P.i_mode < 7
    vnew = t - P.t_stim_sorted(P.i_mode);
    isterminal = 1;
    direction = 1;  %t-tTrigger goes from neg to pos, so 1. 
elseif P.i_mode ==7 %waiting for ground release
    [~,sol]=ode_jumper(t,state,P);
    vnew = sol(5);
    isterminal = 1;
    direction = -1; %ground reaction + to - is -1. 
else
    vnew = t-P.t_max;
    isterminal = 1;
    direction = 1;
end
