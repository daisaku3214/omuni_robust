close all;
clear all;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
time = 2;               %[sec] simuration time
length = time/Deltat.simvel;
robo = omunirobot(parameter,Deltat,[p0;zeros(3+ell,1)],zeros(ell,1));

u0 = [10 -10 -10 10]';
lamp1 = linspace(0,1,length/8);
step = ones(1,6*length/8);
lamp2 = linspace(1,0,length/8);
u = u0*[lamp1 step lamp2];
figure;
for i = 1:(length)/robo.const.ratiovel
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.const.ratiovel
    robo = robo.shiftx;
    robo.u = u(:,robo.const.ratiovel*(i-1)+j);
    end
end
figure;
for i = 1:ell
subplot(ell,1,i);
plot(robo.Tlog,robo.Flog(i,:));
xlim([0 time]);
end
figure;
for i = 1:3
subplot(3,1,i);
plot(robo.Tlog,robo.Xlog(i,:));
xlim([0 time]);
end