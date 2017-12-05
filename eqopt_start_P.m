function [lcerel,gamma,stim,tor,out]=eqopt_start_P(fi,P)
% function [lcerel,gamma,stim,tor,out]=eqopt_start_P(fi,p)
% compute lcerel, gamma, and stim for a starting position in equilibrium,
% given joint angles fi and parameter structure p.
% INPUTS:
% fi rad, external reference frame, from x axis (rightward)
% p parameter structure (getJumperStruct)
% outputs:
% lcerel (dimensionless)
% gamma (dimensionless)
% stim (dimensionless)function [lcerel,gamma,stim,out]=startmP_opt(fi,P)
% 
% NOTE: in startmP_opt we compute torques at the joints if P.tor is empty.

l= P.sk.l(:);
d= P.sk.d(:);
mass=  P.sk.mass(:);
j= P.sk.j(:);

if isfield(P,'tor')
    tor= P.tor;
else
    fprintf('computing g torque since it was not passed to eqopt_start_P.\n');
    tor =get_g_torque(fi,P);
end;
rm= P.m.rm;
c=P.m.c;
q0= P.m.q0;
rldak= P.m.rldak;
ints= P.m.ints;
getal= P.m.getal;
exphat= P.m.exphat;
fmax = P.m.fmax;
rlceopt=  P.m.rlceopt;
width=  P.m.width;
arel=  P.m.arel;
brel=  P.m.brel;
slopfac=  P.m.slopfac;
fasympt=  P.m.fasympt;
sloplin=  P.m.sloplin;
rlsenul=  P.m.rlsenul;
strainmax= P.m.strainmax;
rlpenul=  P.m.rlpenul;
stresmax= P.m.stresmax;
rmapar= P.m.rmapar;

c1= -1.0./(width.^2);
c2= -2.0*c1;
c3=  1.0+c1;
kse=fmax./((strainmax.*rlsenul).^2);
kpe=(stresmax.*fmax)./(((width+1.0-rlpenul).*rlceopt).^2);
kpe = 0;
fprintf('warning! pe has been turned off.\n');

fi(1)=fi(1);
tdiff = diff(fi);
fijo= [fi(1) tdiff(:)'];
%fipjo=[fip(1) diff(fip)];

nmus=length(rlsenul);
nseg = length(P.sk.l);
rloi=zeros(1,nmus);
rmomarm=zeros(nseg,nmus);
for j=1:nseg
    start=(j-1)*3+1;
    rloi=rloi+rmapar(start,:)+rmapar(start+1,:).*fijo(j)+rmapar(start+2,:).*fijo(j)^2;
    rmomarm(j,:)=rmapar(start+1,:)+2*rmapar(start+2,:)*fijo(j);
end


doopt = 1;
if doopt
    opts = optimset('algorithm','sqp','Display','off');
%     opts = optimset('algorithm','active-set','Display','off','GradConstr','On');
    % [fse]=fmincon(@(x)cost_FnormFmax(x,fmax),ones(6,1)*100,[],[],rmomarm,-tor,zeros(6,1),ones(6,1)*10000,[],opts);
    eqopt = struct;
    eqopt.fmax = fmax;
    eqopt.kse =kse;
    eqopt.rlceopt = rlceopt;
    eqopt.rlsenul = rlsenul;
    eqopt.rloi = rloi;
    eqopt.c1 = c1;
    eqopt.c2 = c2;
    eqopt.c3 = c3;
    C_FNORMNOISOM = 1;
    C_FNORM=2;
    C_FNORMLCE = 3; %2017 this is what I've been using. 
    C_FNORMrelFviaLCE = 4;
    %     cmeth=C_FNORM
    cmeth = C_FNORMLCE;
    if cmeth==C_FNORM
        costTxt = 'FNORM';
        ggain = 2;
        [fse,cval,exitflag,output,lambda,grad,hessian]=fmincon(@(x)cost_FnormFmax(x,fmax,eqopt),ones(6,1)*1,[],[],rmomarm,-tor,ones(6,1).*fmax(:).*q0(:)*ggain,ones(6,1).*fmax(:),[],opts);       
    elseif cmeth==C_FNORMNOISOM
        ggain = 10;
        [fse,cval,exitflag,output,lambda,grad,hessian]=fmincon(@(x)cost_FnormFmaxNOISOM(x,fmax,eqopt),ones(6,1)*1,[],[],rmomarm,-tor,ones(6,1).*fmax(:).*q0(:)*ggain,ones(6,1).*fmax(:),[],opts);
        costTxt = 'FNORMNOISOM';
    elseif cmeth==C_FNORMLCE
        x0= ones(6,1)*1000;
        [fse,cval,exitflag,output,lambda,grad,hessian]= ...
            fmincon(@(x)cost_FnormFmaxlce(x,eqopt),...
           x0 ,[],[],rmomarm,-tor(:),zeros(6,1),fmax(:),[],opts);
        [t1,q,cost_f,cost_l,cost_q]=cost_FnormFmaxlce(fse,eqopt);
        fprintf('Stable equilibrium opt: %.2f Force cost= %.2f;Length cost= %.2f; Q cost=%.2f\n',cval,cost_f,cost_l,cost_q);
        costTxt = 'FNORMLCE';
    elseif cmeth == C_FNORMrelFviaLCE
        eqopt.rmomarm = rmomarm;
        eqopt.tor = tor;
        eqopt.fmax = fmax;
        % % % 1- take 1.4 as rellce
        relce_t = ones(6,1)*1.4;
        % % % 2- compute ese as rloi - (rellce)*rlceopt - rlsenul
        ls0s = ones(6,1)*1;
        ls0=[1.15748
            1.05308
            1.32571
            0.948807
            1.08407
            1.25933];
        opts = optimset('algorithm','active-set','Display','off','GradConstr','On');
        
        [rellce,cval,exitflag,output,lambda,grad,hessian]=fmincon(@(x)cost_relFviaLCE(x,fmax,eqopt),...
            [ls0s(:)],[],[],[],[],[ones(6,1)*0.56],[ones(6,1)*1.4],...
            @(x)nlcon_relFviaLCE(x,eqopt),opts);
        cost_relFviaLCE(rellce,fmax,eqopt)
        [a,b,c,d]=nlcon_relFviaLCE(rellce,eqopt)
        rlce = rlceopt(:).*rellce;
        ese = rloi(:)-rlce - rlsenul(:);
        fse = ese.^2.*kse(:);
        %%%we check:
        Fisom_norm = rellce.^2.*eqopt.c1(:)+rellce.*eqopt.c2(:)+eqopt.c3(:);
        q = fse./(Fisom_norm.*fmax(:));
        q(q>1) = 1;
        q(q<0.005) = 0.005;
        fse = fse;
        costTxt = 'FNORMrelFVIALCE';
        
    end;
    % frel = fse./fmax;
    
    fse=fse';
    
else
    torhamknee=-8; %chosen hamstring.
    fse(6)=torhamknee/rmomarm(3,6); %define torque at the knee from hamstring; calculate mforce.
    fse(1)=(-tor(2)+0.5)/rmomarm(2,1); %torque at the ankle is overkill slightly
    fse(2)=-0.5/rmomarm(2,2); %the overkill is precisely matched by the gastrocs. defined gas at knee,rmomarm(3,2)*fse(2).
    fse(3)=(-tor(3)-torhamknee-fse(2)*rmomarm(3,2)-.5)/rmomarm(3,3);
    fse(4)=0.5/rmomarm(3,4);
    fse(5)=(-tor(4)-fse(4)*rmomarm(4,4)-fse(6)*rmomarm(4,6))/rmomarm(4,5);
    costTxt = 'MFB';
    cval = 0;exitflag=0;output=0;lambda=0;grad=0;hessian=0;
    % now we have all the fses, meaning all of the fse lengths.
end;
ese=sqrt((fse)./kse);
ftemp=kse.*ese.^2;
lcerel=(rloi-rlsenul-ese)./rlceopt;
fpeyn=lcerel>rlpenul;
fpe=fpeyn.*(kpe.*((lcerel-rlpenul).*rlceopt).^2);
ftemp=ftemp-fpe;
ftemp=max(ftemp,0*ftemp);

fisom=fmax.*(c1.*lcerel.^2+c2.*lcerel+c3);
q=ftemp./fisom;
% if sum(q>0.99) > 0
%     fprintf('EQ failed to compute active state between 0 and 1.\n');
%     q %show q;
%     stim = zeros(1,6);
%     gamma = zeros(1,6);
% else
%     fprintf('solved starting lengths are %g\n',lcerel);
    q=max(q,ones(size(q)).*q0);
    q=min(q,ones(size(q)).*0.99);
    % % here we catch where there is no stable equilibrium.
%     q=(q0+(rho.*gamma).^exphat)./(1.0+((rho.*gamma).^exphat));
    rho=getal.*((rldak.^ints)-1.0)./(((rldak./lcerel).^ints)-1.0);
    gamma=((q-q0)./(rho.^exphat.*(1-q))).^(1./exphat);
    stim=gamma./c;
    
% end;
out.fse=fse;
out.tor = tor;
out.rmomarm = rmomarm;
out.cost = cval;
out.costType = costTxt;
out.exitflag = exitflag;
out.output = output;
out.lambda = lambda;
out.grad = grad;
out.hessian = hessian;
out.q = q;