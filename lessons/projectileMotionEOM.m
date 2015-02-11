function varargout = projectileMotionEOM(t, x, flag)

    switch flag
        case ''
            varargout{1} = sysEOM(t, x);
        case 'events'
            [varargout{1:3}] = events(x, t(1));
        otherwise
            error(['Unknown flag: ' flag]);
    end
end



function xdot = sysEOM(t, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:     t - vector of times to evaluate the system of equations at
%             x - state space vector
%
% Outputs:    xdot  - solution to set of ODEs
%
% Description: Bouncing ball state space differential equations to be used
%              with ODE solvers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% legend of values in code: 
% 
%       x(1) = horizontal displacement
%       x(2) = vertical displacement
%       x(3) = horizontal velocity
%       x(4) = vertical velocity
%
% same order will be maintained for outputs

xdot = zeros(size(x));
xdot(1) = x(3);
xdot(2) = x(4);
xdot(3) = 0;
xdot(4) = -9.81;

end



function [value,isterminal,direction] = events(x, t0)

value = zeros(1,1);
value(1) = x(2);
isterminal = [1];
direction = [-1];

end