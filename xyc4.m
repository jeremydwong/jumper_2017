function [x,y,xp,yp,xdp,ydp]=xyc4(fi,fip,fidp,xbase,xbasep,xbasedp,l)
% function [x,y,xp,yp,xdp,ydp]=xyc4(fi,fip,fidp,xbase,xbasep,xbasedp,l)
% returns the kinematics of alll bodies' COM. 
% NOTE: expects dimensions for fi,fip,fidp as state x time.
if nargin==2
  l = fip; %this is sloppy. handles passing fi and l for static start use of xyc4. 
  [fip,fidp] = deal(zeros(size(fi,1),size(fi,2)));
  [xbase,xbasep,xbasedp] = deal(zeros(2,size(fi,2)));
end

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
