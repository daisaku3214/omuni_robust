close all;
clear all;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
time = 15;              %[sec] simuration time
times = [1 12 1 1];      %[sec] the time of input shaping
length = time/Deltat.simvel;
robo = omunirobot(parameter,Deltat,[p0;zeros(3+ell,1)],zeros(ell,1));
%input example:ramp input
u0 = v2omega*[-pi/4*2*sin(pi/4) -pi/4*2*sin(pi/4) -pi/4]';
lamp1 = linspace(0,1,times(1)*length);
step1 = ones(1,times(2)*length);
lamp2 = linspace(1,0,times(3)*length);
step2 = zeros(1,times(4)*length);
u = u0*[lamp1 step1 lamp2 step2];


figure;
for i = 1:(length)/robo.const.ratiovel
    %plot trajectory
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    hold on;
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.const.ratiovel
    %plot robot posture
    robo = robo.shiftx;
    robo.u = u(:,robo.const.ratiovel*(i-1)+j);
    end
end
%%
figure;
makemoviefromomunirobot(robo,limit,'movie');
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
