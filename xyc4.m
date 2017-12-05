function [x,y,xp,yp,xdp,ydp]=xyc4(fi,fip,fidp,xbase,xbasep,xbasedp,l)
% function [x,y,xp,yp,xdp,ydp]=xyc4(fi,fip,fidp,xbase,xbasep,xbasedp,l)
% returns the kinematics of the bodies. 
% expects dimensions state x time.

lcfi=diag(l)*cos(fi);
lsfi=diag(l)*sin(fi);

block=...
  [1 0 0 0 0;
   1 1 0 0 0;
   1 1 1 0 0;
   1 1 1 1 0;
   1 1 1 1 1;
  ];
  
x=block*[xbase(1,:);lcfi];
y=block*[xbase(2,:);lsfi];

xp=block*[xbasep(1,:);-lsfi.*fip];
yp=block*[xbasep(2,:); lcfi.*fip];

xdp=block*[xbasedp(1,:);-lsfi.*fidp]+ ...
    block*[xbasedp(1,:)*0;-lcfi.*(fip.^2)];
ydp=block*[xbasedp(2,:); lcfi.*fidp]+ ...
    block*[xbasedp(2,:)*0;-lsfi.*(fip.^2)];
