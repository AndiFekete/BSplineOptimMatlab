function [x,y]=draw_unitcircle()
theta = linspace(0,2*pi,100);
x = cos(theta);                  
y = sin(theta);
plot(x,y, 'LineWidth',4); 
hold on;
plot([-1,1],[0,0],'k',[0,0],[-1,1],'k');
title('Coefficients');
axis('equal'); 
