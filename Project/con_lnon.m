% Non-linear constraints for Non-Linear Constraint Example
function [ c,ceq ] = nonlcon( x )
% Inputs: 
% @x - variable used (i.e. x vector)
% Outputs:
% @c - nonlinear inequality constraints
% @ceq - nonlinear equality constraint function

% Functions must  be in the form f(x) <=0 for c and f(x) = 0 for ceq
% if the Left Side of the equality/inequality is not zero, you must rewrite
% the equation

% Nonlinear constraint for the to do problem

c(1) = @(x) 7*x(1)^3 - 282;
ceq = [];


end
