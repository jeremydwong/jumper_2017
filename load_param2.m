function out = load_param2(param,vmus)
% function out = load_param2(param,vmus)
% out.n and out.p are the parameter structures.
% vmus is the vector of muscles to be used.
% 1   2   3   4    5   6
% MSF MSE MEF MEE  BF  BE

out=struct;
muscles = {'MSF','MSE','MEF','MEE','BF','BE'};

muspar = param(17:50,:);
muspar(35,:) = [1 -1 1 -1 1 -1];
rmapar=   muspar(20:31,:);
segpar = param(1:16,:);

out.n=struct;
out.nmus = length(vmus);
out.n.m = 11.3;%!
out.n.c = 1.37e-4;%!
out.n.eta = 5.27e4;%!%(looks like an 'n'); ACTIVE STATE PARAM
out.n.q0 = 5e-3;%!%active state at time 0? or something?
out.n.k=2.9;%!
out.n.a_rel = 0.41;%!
out.n.b_rel = 5.2;%!
out.n.q_crit = 0.3;
out.n.Ca_gamma_max = 1;%!
out.n.l_pe_f=1.4; %yes this looks similar to l_se_fact; ratio of muscle length
%where PE starts delivering force.s

% location for K_SE measure
out.n.l_se_fact=1.04;%!   %place where k_se is measured.
% parameters for force-velocity
out.n.fasymp = 1.5;%!
out.n.slopfac = 2;%!
out.n.vfactmin = 0.1;%!
out.p = struct;
for i=1:size(vmus,2)
    out.p.muscles{i} = muscles{vmus(i)}; %there's no better way to do this? Ugh?
end;
out.p.FMAX = muspar(8,vmus);%!%1550 N max force of elbow extensor;
out.p.l_ce_opt=muspar(9,vmus);%!
out.p.l_se_0=muspar(16,vmus); %!    %length of tendon at max slacklength, F=0
% the following are l_oi parameters for our elbow flexor.
out.p.a0s=rmapar(7,vmus); %!
out.p.a1s=rmapar(8,vmus); %!
out.p.a2s=rmapar(9,vmus);%!
out.p.a0e=rmapar(10,vmus);%!
out.p.a1e=rmapar(11,vmus);%!
out.p.a2e=rmapar(12,vmus);%!
out.p.width = muspar(10,vmus);%! %operating width of muscles about normalized optim length of 1
%                 this is actually not the case for shoulder muscles, which
%                 conform to the sliding-filament-theory .56 (since they
%                 are not lumped muscles)
out.p.a=1./(out.p.width).^2; %for calculating force-length relationship.
out.sk.g = 9.81*0;
% out.sk.m1 = 2.5;
% out.sk.m2 = 2.5;
% out.sk.l1 = 0.3;
% out.sk.l2 = 0.3;
% out.sk.I1 = (1/12)*out.sk.m1*out.sk.l1*out.sk.l1;
% out.sk.I2 = (1/12)*out.sk.m2*out.sk.l2*out.sk.l2;
% out.sk.r1=out.sk.l1*0.5;
% out.sk.r2=out.sk.l2*0.5;
out.sk.m1 = 2.1;
out.sk.m2 = 1.65;
out.sk.l1 = .3348;
out.sk.l2 = 0.2628;
out.sk.r1 = 0.1460;
out.sk.r2 = 0.1792;
out.sk.I1 = 0.0244;
out.sk.I2 = 0.0250;

out.nmus = 6;
out.njoints = 2;

%keep the following to matchw with the old code.
% out.n.m1=2.5; %random guess at lower arm mass; 3% of total body weight 70 kilo.
% out.n.l1 = 0.3;
% out.n.I1 = (1/12)*out.n.m1*out.n.l1*out.n.l1;


% out.p.dir=muspar(35,vmus);
% EE    
% out.p.a0=0.236;
% out.p.a1e=0.025;
% out.p.a1s=0.0;
% out.p.a2e=-2.16e-3;
