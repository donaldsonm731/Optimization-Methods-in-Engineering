%% Problem 3

f_star = 1000;

for d =  0.015:0.015:0.25
    for D = 1/8:1/8:1.5
        for N = 2:100

            x = [d,D,N];
            [f,g1,g2,g3,g4] = p3(x);
            if g1 < 0 && g2 < 0 && g3 < 0 && g4 < 0 && f < f_star
                f_star = f;
                x_star = x;
            end

        end
    end
end

f_star
x_star
%% Problem 4
clc
clear

% Code from mathworks:
% https://www.mathworks.com/matlabcentral/answers/400426-optimization-of-objective-function-with-multiple-constraints-and-variables

fun=@(x) (x(1)-10)^2 + 5*(x(2) - 12)^2 + 3*(x(4) - 11)^2 + 10*x(5)^6 ...
          + 7*x(6)^2 + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);

c = @(x) [2*x(1)^2 + 3*x(2)^4 + x(3) + 4*x(4)^2 + 5*x(5) - 127;

7*x(1) + 3*x(2) + 10*x(3)^2 + x(4) - x(5) - 282;

23*x(1) + x(2)^2 + 6*x(6)^2 - 8*x(7) - 196;

4*x(1)^2 + x(2)^2 - 3*x(1)*x(2) + 2*x(3)^2 + 5*x(6) - 11*x(7);

]; 

ceq = @(x)[];
nonl_con = @(x)deal(c(x),ceq(x));

gs = GlobalSearch;
opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','x0',[2,2,1,4,0,1,2],'objective',fun,'lb',[2,2,2,-Inf,-Inf,-Inf,-Inf],'ub',[2,2,2,Inf,Inf,Inf,Inf],'nonlcon',nonl_con,'options',opts);

min=run(gs,problem);
min
x = min;
fun(x)
c(x)

% %x = [5,5,5,11,0,1.8,3.8];
% %f_star = 1000;
% %x0 = [0,0,0,0,0,0,0];
% 
% f = @(x) (x(1)-10)^2 + 5*(x(2) - 12)^2 + 3*(x(4) - 11)^2 + 10*x(5)^6 ...
%           + 7*x(6)^2 + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);
% % x0 = [5,5,5,10,10,10,10];
% % A = [];
% % b = [];
% % Aeq = [];
% % beq = [];
% % lb = [-Inf,-Inf,-Inf, -Inf, -Inf, -Inf, -Inf];
% % ub = [Inf,Inf,Inf, Inf, Inf, Inf, Inf];
%  %nonlcon = [];
%  %nonlcon = nonlcon;
%  %x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,'nonlcon');
% 
% 
% [x,fval,exitflag,output] = fminsearch(f,x0)
% 
% g1 = @(x) 2*x(1)^2 + 3*x(2)^4 + x(3) + 4*x(4)^2 + 5*x(5) - 127;
% g2 = @(x) 7*x(1) + 3*x(2) + 10*x(3)^2 + x(4) - x(5) - 282;
% g3 = @(x) 23*x(1) + x(2)^2 + 6*x(6)^2 - 8*x(7) - 196;
% g4 = @(x) 4*x(1)^2 + x(2)^2 - 3*x(1)*x(2) + 2*x(3)^2 + 5*x(6) - 11*x(7);
% f_star = 10000000000000000000;
% for i = 1:5
%     for j = 1:5
%         for k = 1:5
%             for l = 0:10
%                 for m = 0:10
%                     for n = 0:10
%                         for o = 0:10
% 
%                            x = [i,j,k,l,m,n,o];
% 
%                            if f(x) < f_star && g1(x) < 0 && g2(x) < 0 && g3(x) < 0 && g4(x) < 0
%                                 f(x)
%                                 x
%                                 f_star = f(x);
%                                 x_star = x;
%                            end
%                            
% 
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end



%% Problem 5
clc
clear all

w = [19 15 20 8 5 7 3 2 4];
v = [380 225 320 96 70 126 30 22 68];
value_star = 0;

for i = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                for m = 0:1
                    for n = 0:1
                        for o = 0:1
                            for p = 0:1
                                for q = 0:1
                                    weight = sum(w.*[i,j,k,l,m,n,o,p,q]);
                                    value = sum(v.*[i,j,k,l,m,n,o,p,q]);

                                    if weight <= 40 && value > value_star
                                        value_star = value;
                                        x_star = [i,j,k,l,m,n,o,p,q];
                                        weight_star = weight;

                                    end
                                     

                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

weight_star
x_star



%% Problem 6
clc
clear all

w = [20, 12, 7, 75, 93, 21, 75, 67, 34, 28];
vol = [41, 51, 24, 40, 84, 70, 34, 41, 49, 27];
val = [84, 34, 31, 14, 67, 65, 86, 98, 50, 7]; 
value_star = 0;

for i = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                for m = 0:1
                    for n = 0:1
                        for o = 0:1
                            for p = 0:1
                                for q = 0:1
                                    for r = 0:1
                                        knapsack = [i,j,k,l,m,n,o,p,q,r];
                                        weight = sum(w.*knapsack);
                                        volume = sum(vol.*knapsack);
                                        value = sum(val.*knapsack);

                                        if weight <= 190 && volume <= 250 &&value > value_star
                                            value_star = value;
                                            x_star = [i,j,k,l,m,n,o,p,q];
                                            weight_star = weight;
                                            volume_star = volume;

                                        end
                                     
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

x_star
weight_star
volume_star







%%



%%
x = min;
fun(x)
c(x)

function [f,g1,g2,g3,g4] = p3(x)
d = x(1);
D = x(2);
N = x(3);

% constants
Q = 2;
P = 10;
g = 386;
rho = 7.38342*10^(-4);
G = 1.15*10^7;

t_a = 80000;
delta = 0.5;
w_0 = 100;
D_0 = 1.5;

K = (d^4*G)/(8*D^3*N);

k = (4*D - d)/(4*(D-d))+ 0.615*d/D;
t = (8*k*P*D)/(pi*d^3);

w = d/(2*pi*N*D^2)*sqrt(G/(2*rho));


f = (N + Q)*D.*d.^2;

g1 = -P/K + delta; % Deflection limit
g2 = t - t_a; % Shear Stress
g3 = w_0 - w ; % Freq. surge waves
g4 = D + d - D_0; % Outer Diameter of spring

end

