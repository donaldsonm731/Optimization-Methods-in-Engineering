%%% MEMS 5001: Hw 5 %%%
% Matthew Donaldson
% April 4, 2022
%% Problem 1
% Part A: GA
% Obj. Function
f = @(x) 100*(x(1) - x(2))^2 + (1-x(1))^2;

% Linear Inequality Constraints
A = [];
b = [];

% Linear Equality Constraints
Aeq = []; 
beq = [];

% Boundry Conditions
LB = [-2.048 -2.048];
UB = [2.048 2.048];
nonlcon = [];

nvars = 2;
opts = optimoptions(@ga,'PlotFcn', {@gaplotbestf});
[x_ga,fval_ga,exitflag_ga,output_ga,population_ga] = ga(f,nvars,A,b,Aeq,beq,LB,UB,nonlcon,opts);

x_ga
fval_ga
% Part B: Design Space
[x1,x2] = meshgrid(-2.048:0.01:2.048,-2.048:0.01:2.048);
func = 100*(x1-x2).^2 + (1-x1).^2;

figure(1);
contour(x1,x2,func, 20);

figure (2);
mesh(x1,x2,func)

%% Probem 2

% Part A: GA
% Obj. Function
f2 = @(x) 2.*x(1).^3 + 15.*x(2).^2 - 8.*x(1).*x(2) - 4.*x(1);

% Linear Inequality Constraints
A = [];
b = [];

% Linear Equality Constraints
Aeq = []; 
beq = [];

% Boundry Conditions
LB = [-4 -8];
UB = [4 8];
nonlcon = @hw5NoLinearIneq;

nvars = 2;
opts = optimoptions(@ga,'PlotFcn', {@gaplotbestf});
[x2_ga,fval2_ga,exitflag_ga,output_ga,population_ga] = ga(f2,nvars,A,b,Aeq,beq,LB,UB,nonlcon,opts);

x2_ga
fval2_ga


% Part B: Design Space
[x1,x2] = meshgrid(-4:0.01:4,-8:0.01:8);
func2 = 2*x1.^3 + 15*x2.^2 - 8*x1.*x2 - 4*x1;

figure(1)
contour(x1,x2,func2, 20);
figure(2)
mesh(x1,x2,func2)


%% Problem 3

% Part 1
% Obj. Function
 f = @(x) (-1*(1*(x-2).^2 + 2) + 45)*(0 >= x & x< 4)...
        + (-1*(10*(x-6).^2 + 1) + 45)*(4>= x & x < 8)... 
        + (-1*(2*(x-12).^2 + 10) + 45)*(8>= x & x < 16)...
        + (-1*(0.1*(x-24).^2 + 1.1) + 45)*(16>= x & x <= 31);

% Part 2
fplot(f)

% Part 3
% Linear Inequality Constraints
A = [];
b = [];

% Linear Equality Constraints
Aeq = []; 
beq = [];

% Boundry Conditions
LB = [0];
UB = [31];
nonlcon = [];

nvars = 1;
opts = optimoptions(@ga,'PlotFcn', {@gaplotbestf});
[x_ga,fval_ga,exitflag_ga,output_ga,population_ga] = ga(f,nvars,A,b,Aeq,beq,LB,UB,nonlcon, opts);

x_ga
fval_ga



%x = 0:0.01:31;
%func = (-1*(1*(x-2).^2 + 2) + 45).*(0 >= x & x< 4)...
%         + (-1*(10*(x-6).^2 + 1) + 45).*(4>= x & x < 8)... 
%         + (-1*(2*(x-12).^2 + 10) + 45).*(8>= x & x < 16)...
%         + (-1*(0.1*(x-24).^2 + 1.1) + 45).*(16>= x & x <= 31);



% Part of Problem 2: Nonlinear constraint
% Non-linear constraints for Non-Linear Constraint Example
function [ c,ceq ] = hw5NoLinearIneq( x )
% Inputs: 
% @x - variable used (i.e. x vector)
% Outputs:
% @c - nonlinear inequality constraints
% @ceq - nonlinear equality constraint function

% Functions must  be in the form f(x) <=0 for c and f(x) = 0 for ceq
% if the Left Side of the equality/inequality is not zero, you must rewrite
% the equation

% Nonlinear constraint for the to do problem
c(1) = x(1) - 0.25.*x(2).^2 - 1;
c(2) = -x(1) - x(2).^2;

ceq = [];


end







