
function animateJump(x,y,cmx,cmy,fname)
figure;
ds = 30;
disp('press any key');
%downsample
i_downsamp=1:ds:length(x);
x = x(i_downsamp,:);
y = y(i_downsamp,:);
cmx = cmx(i_downsamp);
cmy = cmy(i_downsamp);
ax = [-1,1,0,2];
F(length(cmx)) = struct('cdata',[],'colormap',[]);
for ib=1:length(cmx)
    hold off;
    
    plot(x(ib,:),y(ib,:),'-');
    hold on;
    plot(cmx(ib),cmy(ib),'bo','markersize',20);
        
    hold off;
    axis(ax);
    drawnow;
    pause(1/ds);
    F(ib) = getframe(gcf);
end;

if nargin ==5
    v = VideoWriter(fname,'MPEG-4');
    %v.CompressionRatio=10;
%     v.FileFormat='mp4';
    open(v);
    writeVideo(v,F);
    close(v);
end;