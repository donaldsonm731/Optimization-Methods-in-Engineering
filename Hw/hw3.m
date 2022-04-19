% Problem 6

[r, h] = meshgrid(0:0.01:20, 0:0.01:20);

f = pi.*r.^2.*h;
g1 = -2*pi.*r.*h - 900;
g2 = r - 20;
g3 = -r - 5;
g4 = h - 20;
g5 = h;


f_bar = -540*pi - 90*pi*(r-6) - 36*pi*(h-15);
g1_bar = -180*pi - 900 - 30*pi*(r-6) - 12*pi*(h-15);
g2_bar = g2;
g3_bar = g3;
g4_bar = g4;
g5_bar = g5;

contour(r,h,f,30)
hold on 
plot(6,15,'*k','linewidth', 3)
hold on
contour(r,h,g1,[0 0],'-k', 'LineWidth', 3);


%contour(r,h,g2,[0 0], 'k','LineWidth', 3)
%hold on 
%contour(r,h,g3,[0 0], 'k','LineWidth', 3)
%hold on 
%contour(r,h,g4,[0 0], 'k','LineWidth', 3)
%hold on 
%contour(r,h,g5,[0 0], 'k','LineWidth', 3)
hold on
contour(r,h,f_bar,10,'--r');
hold on
contour(r,h,g1_bar,[0 0],'--r', 'LineWidth', 3)

title('Matthew Donaldson')
xlabel('r')
ylabel('h')
axis([5 20 0 20])
legend('Original f','Point of Interst', 'Surface Area Constraint', 'Linearized f', 'Finearized Surface Area Constraint')