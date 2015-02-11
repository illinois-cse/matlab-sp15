function dN = decay_scheme(t,N0)
% We expect that N0 is a six-element vector.

% Half-lives / s
t_Rn = 5.4e4;
t_At = 2.6e4;
t_Po211 = 0.52;
t_Po207 = 2e4;
t_Bi = 9.46e8;

% Decay constants / s^-1
l_Rn = log(2)/t_Rn;
l_At = log(2)/t_At;
l_Po211 = log(2)/t_Po211;
l_Po207 = log(2)/t_Po207;
l_Bi = log(2)/t_Bi;

% Equation matrix
eqns = ...
    [-l_Rn 0 0 0 0 0;
     0.74*l_Rn -l_At 0 0 0 0;
     0 l_At -l_Po211 0 0 0;
     0.26*l_Rn 0 0 -l_Po207 0 0;
     0 0 0 l_Po207 -l_Bi 0;
     0 0 l_Po211 0 l_Bi 0];
 
 % Differential change from N0
 dN = eqns*N0;