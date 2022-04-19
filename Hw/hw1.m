%% Beer can Problem
clf
clc
% Mesh grid 
[D,H] = meshgrid(3.5:0.01:8,8:0.01:18 );

% Obj. Function we are trying to minimize
A = pi*D.*H + pi*D.^2/4;

% Constraints
hMax = 18; %g1
hMin = 8; % g2
dMax = 8; % g3
dMin = 3.5; % g4
h = 1600./(pi*D(1,:).^2); %g5

% Plotting the contour with all the boundry conditions
[C,z] = contour(D,H,A,20);
hold on 
yline(hMax, 'k', 'linewidth',3)
yline(hMin, 'b', 'linewidth',3)
xline(dMax, 'g', 'linewidth',3)
xline(dMin, 'r', 'linewidth',3)
hold on
plot(D(1,:),h, 'k--', 'linewidth', 2)
hold on
plot(8, 8, 'ro', 'LineWidth', 2, 'MarkerSize', 12);

% Labels and axis limits
text(6.2,17.5,'Feasable Region')
% text(6.6,18.3,'g1')
% text(6.6,7.8,'g2')
% text(8.08,13,'g3')
% text(3.1,12,'g4')
% text(6.15,13.75,'g5')
% str = {'\leftarrow Optimal','       Value'};
% text(8.1,8,str)

clabel(C,z)
axis([5 9 7.5 18.5])
title('Matthew Donaldson')
xlabel('Diameter')
ylabel('Height')
legend('Contour Lines','g1', 'g2', 'g3', 'g4','g5','Optimal Value' )
