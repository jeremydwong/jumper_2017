function [cmx,cmy,cmxp,cmyp,cmxdp,cmydp]=cm4(x,y,xp,yp,xdp,ydp,l,d,m)
% function [cmx,cmy,cmxp,cmyp,cmxdp,cmydp]=cm4(x,y,xp,yp,xdp,ydp,l,d,m)
% compute kinematics of centre of masses of all segments.

if(l(4)==0)
  l(4)=1e-10;
  d(4)=1e-10;
end
block= [   1     d(1)/l(1)    0       0        0        0        0        0    ;
           0         0        1   d(2)/l(2)    0        0        0        0    ;
           0         0        0       0        1   d(3)/l(3)     0        0    ;
           0         0        0       0        0        0        1   d(4)/l(4)];
         
         
cmx=[m/sum(m)]'*block*[x(1,:);x(2,:)-x(1,:);x(2,:);x(3,:)-x(2,:);x(3,:);x(4,:)-x(3,:);x(4,:);x(5,:)-x(4,:)];
cmy=[m/sum(m)]'*block*[y(1,:);y(2,:)-y(1,:);y(2,:);y(3,:)-y(2,:);y(3,:);y(4,:)-y(3,:);y(4,:);y(5,:)-y(4,:)];
cmxp=[m/sum(m)]'*block*[xp(1,:);xp(2,:)-xp(1,:);xp(2,:);xp(3,:)-xp(2,:);xp(3,:);xp(4,:)-xp(3,:);xp(4,:);xp(5,:)-xp(4,:)];
cmyp=[m/sum(m)]'*block*[yp(1,:);yp(2,:)-yp(1,:);yp(2,:);yp(3,:)-yp(2,:);yp(3,:);yp(4,:)-yp(3,:);yp(4,:);yp(5,:)-yp(4,:)];
cmxdp=[m/sum(m)]'*block*[xdp(1,:);xdp(2,:)-xdp(1,:);xdp(2,:);xdp(3,:)-xdp(2,:);xdp(3,:);xdp(4,:)-xdp(3,:);xdp(4,:);xdp(5,:)-xdp(4,:)];
cmydp=[m/sum(m)]'*block*[ydp(1,:);ydp(2,:)-ydp(1,:);ydp(2,:);ydp(3,:)-ydp(2,:);ydp(3,:);ydp(4,:)-ydp(3,:);ydp(4,:);ydp(5,:)-ydp(4,:)];

% cmx = 
