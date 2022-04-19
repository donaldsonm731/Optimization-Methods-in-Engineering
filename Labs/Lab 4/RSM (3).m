% ********************************************
% Suppose f is a function of two design variables x1 and x2 as follows:
%             f = 3*x1^3+10*x2^2-4*x1*x2-6*x1
%
% Generate a response surface with all of the linear and quadratic terms 
% using the least squares method. (number of sample points = 9)
% ********************************************

clc
clear

% We want to approximate evaluations of our problem via the form:
% f_hat = d1 + d[2:n]*z
% where d is a coefficient vector and z is in this case, the linear and
% quadratic terms of x
%
% We solve for the d coefficients though the linear system of equations:
% A*d = b
% We Initialize A and b as empty vectors/matrices of proper size
% For a two variable linear and quadratic term approximation, we will have
% 6 terms in the b vector, and a 6x6 A matrix (yielding a vector of size 6
% for d)
A=zeros(6);
b=zeros(6,1);
% k: number of sample points
k=9;
% sample points X=(x1,x2)
x1=[-1.5,-1.5,-1.5,1.25,1.25,1.25,4,4,4];
x2=[-3,0,3,-3,0,3,-3,0,3];
X=[x1;x2];

% Define the linear and quadratic terms based on sampled points (x1, x2)
% z1 = x1, z2 = x2, z3 = x1^2, z4 = x2^2, z5 = x1*x2
z1=x1;
z2=x2;
z3=x1.^2;
z4=x2.^2;
z5=x1.*x2;
z=[z1;z2;z3;z4;z5]'
% Show all Linear and Quadratic Terms
z

%Calulate A
A(1,1)=k;
A(1,2:6)=sum(z);
A(2:6,1)=sum(z)';
for i=2:6
    for j=2:6
        A(i,j)=sum(z(:,i-1).*z(:,j-1));
    end
end
%Show A matrix
A

% Calculate b
f = @(x) 3*x(1)^3+10*x(2)^2-4*x(1)*x(2)-6*x(1);
for i=1:k
    F(i)=f(X(:,i));
end
b(1)=sum(F);
for i=1:5
    b(i+1)=sum(F'.*z(:,i));
end
% Show b coefficients
b
% solve linear equation A*d=b
d = linsolve(A,b);
% Show d coefficients
d
%Response Surface Model is an approximation for f that uses d as the
%coefficients for the linear and quadratic order terms (z vector) above.
% You could use a cubic or higher order RSM could use higher order terms as well
% NOTE: d1 is an intercept, not a coefficients

% Compare approximation for sample points with true values
approx = d(1) +z* d(2:6)
table(F',approx)
% With a limited number of points (such as in this case), the approximation
% should very closely or exactly fit the true values. However, in general
% we are using a least squares approximation, so this is not generally:

% Try varying x1 and x2 above, and then comparing the true and approximate
% values again

% Try with test points
x1_new = 8
x2_new = 8
z_test = [x1_new x2_new x1_new.^2 x2_new.^2 x1_new.*x2_new]
approx_check = d(1) +z_test*d(2:6)
% Evaluate the true function at this point. Is it the same

