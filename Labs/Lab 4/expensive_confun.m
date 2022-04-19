function [c,ceq] = expensive_confun(x)
%EXPENSIVE_CONFUN An expensive constraint function used in optimparfor example.

%   Copyright 2007-2013 The MathWorks, Inc.

% Simulate an expensive function by pausing
pause(0.1);
% Evaluate constraints
c = [x(1)*x(2)*x(3) - x(1) - x(2) - -x(3) - x(4) 
      + 8];
% No nonlinear equality constraints:
ceq = [];