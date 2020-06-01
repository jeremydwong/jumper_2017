% produce single jump from starting initial condition.
% clear all;clc;
plot_lines = {'linewidth',2};
col=get(groot,'DefaultAxesColorOrder');
%
P = getJumperParametersJDW();


%initial conditions
fi0=[2.5277    0.8295    2.5385    0.7504];

GO_OPTIMIZE_STIM = 1;
if GO_OPTIMIZE_STIM
  timeStimFull =[];
  heightBest = 0;
  nTries=2;
  for i =1:nTries
    x0=ones(1,6)*.04+rand(1,6)*0.01;
    [currentTiming,currentHeight] = fminsearch(@(x)run4seg_P_optstart(fi0,x,P,0),x0);
    if currentHeight < heightBest
      timeStimFull = currentTiming;
      heightBest = currentHeight;
    end
    heightsAll(i) = currentHeight;
    timesAll(i,:) = currentTiming(:)';    
  end
else %use some pre-computed onset times that are good!
  
  timeStimFull= [  0.0716
    0.0746
    0.0182
    0.0462
    0.0030
    0.0487];
end
%% run the solution!

[height,state,fwdData]=run4seg_P_optstart(fi0,...
timeStimFull,P,0);

e = energy(state,fwdData,P);

animate_jump(fwdData.x,fwdData.y,fwdData.cmx,fwdData.cmy);