function [x,y,dx,dy,ddx,ddy]=centreOfMass(fi,P)
% function [x,y,dx,dy,ddx,ddy]=centreOfMass(fi,P)
% compute COM for given fi and param structure P. 
% NOTE: fi should be state x time (length(time) = 1 for use in equilibrium)

[x4,y4,dx4,dy4,ddx4,ddy4]=xyc4(fi,P.sk.l);
[x,y,dx,dy,ddx,ddy]=xyCOM(x4,y4,dx4,dy4,ddx4,ddy4,P.sk.l,P.sk.d,P.sk.mass);
