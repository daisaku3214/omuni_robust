clear all;
close all;
robot_define;
robo = omunirobot(parameter,Deltat,zeros(6+ell,1),zeros(ell,1));
R = 3;
Omega = 1;
T = 2*pi/Omega;
limit = [-5 5; -5 5];
numt = T/Deltat.simpos;
T = linspace(0,T,numt);
X = R.*cos(Omega.*T);
Y = R.*sin(Omega.*T);
Theta = Omega.*T;
xx = [X;Y;Theta];
figure;
for  i = 1:length(T)
 plot(X(:,1:i),Y(:,1:i),'b--');
 hold on;
 robo.x(1:3) = xx(:,i);
 plot_omunirobot(robo,limit);
 hold off;
 drawnow;
end