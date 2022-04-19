objFun =@(x) (x(1)-10)^2 + 5*(x(2) - 12)^2 + 3*(x(4) - 11)^2 + 10*x(5)^6 ...
          + 7*x(6)^2 + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);

% Define an inital guess for x
x0 = [0 0];

% Define the rest of the input arguements as empty vectors
A = [];
b = [];
Aeq = [];
beq = [];
LB = [];
UB = [];
nonlcon = [];

% Call the nonlinear program solver
[x, FunValue, ExitFlag] = fmincon(objFun, x0, A,b,Aeq,beq,LB,UB,nonlcon)