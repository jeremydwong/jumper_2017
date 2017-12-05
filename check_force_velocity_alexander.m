% check_force_velocity_alexander
% lengthening
% shortening is negative length change
v_max = 1;
Fm_iso_all = [1];
% v = -v_max:.01:1.8;
v = -v_max:.01:0
for i=1:length(Fm_iso_all)
    Fm_iso = Fm_iso_all(i);
    Fm = (-v < 0) .* (Fm_iso * (1.8 - .8*(v_max - v)./(v_max + 23*v))) ... %eccentric arm
        + (-v >= 0) .* (Fm_iso * (v_max + v)./(v_max - 3*v)); %concentric arm.
    
    
    % inverses are
    k0 = Fm/Fm_iso - 1.8;k1 = Fm/Fm_iso;
    v_inv = (-v < 0) .* ((k0+0.8)*v_max ./ (0.8 - k0*23)) ...
        + (-v >= 0) .* ((k1-1)*v_max./(1+3*k1));
    plot(v_inv,Fm,'linewidth',2); hold on;
end
grid on;
 axis([-1,1,0,2]);
 figurefyTalk();
xlabel('normalized velocity (v/v_{max})','interpreter','tex');
ylabel('normalized force')


