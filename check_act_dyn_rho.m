m=11.3;c = 1.37e-4;n=5.27e4;
k = 2.9;

lce = 0.66:0.01:1.44;

rho = c*n*lce.*(k-1)./(k-lce)
plot(lce,rho)
xlabel('lce-norm');
ylabel('rho');
%%
rhos = [4;14];
rho = rhos(2);
q0 = 0.005;
gamma = 0.005:0.01:1;
exphat=3;
q = (q0+(rho*gamma).^exphat)./(1+(rho*gamma).^exphat);
hold on;plot(gamma,q);xlabel('gamma ([Ca^{2+}])');ylabel('q (active state)');