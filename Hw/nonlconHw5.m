% Non-linear constraints for Non-Linear Constraint Example
function [ c,ceq ] = nonlconHw5( x )
% Inputs: 
% @x - variable used (i.e. x vector)
% Outputs:
% @c - nonlinear inequality constraints
% @ceq - nonlinear equality constraint function

% Functions must  be in the form f(x) <=0 for c and f(x) = 0 for ceq
% if the Left Side of the equality/inequality is not zero, you must rewrite
% the equation

% Nonlinear constraint for the HW5 P2
c(1) = @(x) x(1) - 1/4*x(2)^2 - 1;
c(2) = @(x) -x(1) - x(2)^2;

ceq = [];


end
