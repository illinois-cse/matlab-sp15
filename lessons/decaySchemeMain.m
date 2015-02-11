% THIS TAKES FOREVER... better to use ode15s
% tic
% [t, N] = ode45(@decay_scheme, [0 10e10], [100;0;0;0;0;0]);
% toc

tic
[t,N] = ode15s(@decay_scheme, [0 10e10], [100;0;0;0;0;0]);
toc

semilogx(t,N)
legend('Rn_211','At_211','Po_211','Po_207','Bi_207','Pb_207')