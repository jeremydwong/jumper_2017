function out = overwriteparamsRLC201306(b)
% function out = overwriteparamsRLC201306(b,newb)
% miraculously, rmapar does not change. 
out = b;
fprintf('sk.l:');
out.sk.l
% L: JumperParams is 0.1650    0.4580    0.4851    0.2922
out.sk.l = [0.1650,0.4449,0.4451,0.5];
out.sk.l
fprintf('sk.d:');
out.sk.d
out.sk.d = [.12,.2523,.2524,.2922];
% out.sk.d = [.12,.2523,.2524,.292168];
% d: JumperParams is 0.1200    0.2597    0.2750    0.2922
out.sk.d
fprintf('m.getal:');
out.m.getal
out.m.getal = repmat(66200,1,6); %%from 52700
out.m.getal
fprintf('m.rlpenul:');
out.m.rlpenul
out.m.rlpenul = repmat(1.5,1,6); %%from 1.4, for opts.
out.m.rlpenul

fprintf('m.exphat:');
out.m.exphat
out.m.exphat = repmat(2,1,6); %%from 3
out.m.exphat

% out.sk.kjo = [0;700;1000;1000];
% out.sk.djo = [0;20;40;40];

out.sk.fiRange=[pi/2,2.7438;
        -2.2,-0.1;
        .2,2.7;
        -3,-.1];
%rename
% out.fmax = gmax;