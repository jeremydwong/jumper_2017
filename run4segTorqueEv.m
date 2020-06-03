function [height,state,fwd]=run4segTorqueEv(in_state,P,varargin)

% tor_spline = [];
% for i_tor = 1:size(length(a'))
%     tor_
% end;
% tor = ppval(t_and_stim(:,1),t_and_stim(:,2:end));
% t_cur =0;

%%TEMPLATE VARARGIN
%%only 1 optional argument, doFlight.
numvarargs = length(varargin);
num_maxvarargs = 1;
if numvarargs > num_maxvarargs
    error(['scoreETF:TooManyInputs', ...
        'requires at most ',num2str(num_maxvarargs),' optional inputs']);
end
% set defaults for optional inputs
doFlight0=0;
optargs = {doFlight0};
% now put these defaults into the optargs cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names used below.
[doFlight] = optargs{:};
% %%/TEMPLATE VARARGIN

air = 0;
P.sim.air = air;

fwd = struct;
fwd.state = in_state(:)';
fwd.state_dt = zeros(1,12);
fwd.t = 0;

for iu =1:length(P.U)
    fwd.tor(iu) = ppval(P.U{iu},0);
end

% events
P.i_mode = 1;
P.DO_FLIGHT = 0;
P.t_max = .7;
odeopts = odeset('events',@eventsTorqueJumper);
% /events

ode_handle = @odeTorqueJumper;P.ode_handle = ode_handle;
solver_handle = @ode45;
ct = 0; stepsize = 0.001;
state0 = in_state;
t_all = [];
state_all = [];
while(P.i_mode < 2 & ct<P.t_max) % loop through timesteps and simulate.
    [t,state,te,ye,ie]=solver_handle(ode_handle,ct:stepsize:P.t_max,state0,odeopts,P);
    t_all=[t_all;t(1:end-1)];
    state_all = [state_all;state(1:end-1,:)];
    % note: annoyingly, removing the end doesn't guarantee you don't get
    % duplicate times. this can still happen it appears in the case where
    % a triggered event happens exactly on your return time index
    % frequency, or in this case, every 1000th of a second.
    
    ct = t(end);
    state0 = state(end,:);
    P.i_mode = P.i_mode + 1; %shift mode schedule enum up 1.
end
t_all = [t_all;t(end)];
state_all = [state_all;state(end,:)];
state = state_all;
t = t_all;
nseg =4;

%% compute jump height.
fi=      state(:,1:nseg);
fip=     state(:,nseg+1:2*nseg);
xbase=   state(:,2*nseg+1:2*nseg+2);
xbasep=  state(:,2*nseg+3:2*nseg+4);
[x,y,xp,yp,~,~]=xyc4(fi',fip',zeros(size(fip))',xbase',zeros(size(xbase))',zeros(size(xbase))',P.sk.l);
[cmx,cmy,cmxp,cmyp,~,~]=kinematics4com(x,y,xp,yp,zeros(size(xp)),zeros(size(yp)),P.sk.l,P.sk.d,P.sk.mass(:));

height=cmy(end)+0.5/9.81*cmyp(end)^2;
height = -height;
%% /compute jump height.

if nargout == 3
    [sol,tor]=deal([]);
    for ifwd =1:length(t_all)
        [stated_cur,sol_cur,os]=ode_handle(t(ifwd),state(ifwd,:),P);
        sol = [sol; sol_cur'];
        tor = [tor; os.tor(:)'];
    end;
    fwd.sol = sol;
    fwd.tor = tor;
    fwd.phi = state_all(:,1:4);
    fwd.phidot = state_all(:,5:8);
    fwd.t = t_all;
    fwd.x = x';
    fwd.y = y';
    fwd.cmx = cmx';
    fwd.cmy = cmy';
end;
%catch too long muscles
fprintf('height = %.3f \n',height);
