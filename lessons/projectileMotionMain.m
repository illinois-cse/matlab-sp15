%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Erika Fotsch
%  Date:   1/22/2015
%
%  Description:  this example will show how to set up event handling with
%  Matlab's built-in ode solvers.  This file will be the script to run and
%  another file will be created which will house the differential equaitons
%  to solve and the event crossing function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace and variables
clear;
clc;

% set initial conditions
t0 = 0;                 % initial time = 0 sec.
v0 = 0;                 % initial velocity
theta = 30*pi/180;      % angle of projection
x = zeros(1,4);
x(1) = 0;
x(2) = 1000;
x(3) = v0*cos(theta);
x(4) = v0*sin(theta);
tf = 100;               % simulation run time = 100 sec.
e = 0.8;                % coefficient of restitution for impact law


% initialize arrays to store all data
X = [];
T = [];


% turn on events for ode solver
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'events', 'on');


while t0 < tf
    
    [tout, xout, te, xe, ie] = ode45('projectileMotionEOM', [t0 tf], x, options);
    
    X = [X; xout];
    T = [T; tout];
    
    if tout(end) == tf
        break;
    end
    
    % set new initial conditions
    x = xout(end,:);
    t0 = tout(end);
    
    if ~isempty(ie)
       if ie(end) == 1
           x(4) = -e*x(4);
           stat = 1;
       end
    end
end

plot(T,X(:,2));