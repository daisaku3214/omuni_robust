close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
robo = omunirobot(parameter,Deltat,Matrix,[p0;zeros(3+ell,1)],zeros(ell,1));
%input example:ramp input:move circle
% u1 = v2omega*[-pi/4*2 0 -pi/4]';
% times = [2 10 2 1];     %[sec] the time of input shaping
% time = sum(times);      %[sec] simuration time
% length = time/Deltat.simvel;
% lamp1 = linspace(0,1,times(1)/time*length);
% step1 = ones(1,times(2)/time*length);
% lamp2 = linspace(1,0,times(3)/time*length);
% step2 = zeros(1,times(4)/time*length);
% u = u1*[lamp1 step1 lamp2 step2];
%input example:ramp input:move linear
times = [1 1 1 1 1 1 1 1];     %[sec] the time of input shaping
time = sum(times);      %[sec] simuration time
length = time/Deltat.simvel;
u1 = v2omega*[2 0 0]';
u2 = v2omega*[-2 -2 0]';
lamp1 = linspace(0,1,times(1)/time*length);
step1 = ones(1,times(2)/time*length);
lamp2 = linspace(1,0,times(3)/time*length);
step2 = zeros(1,times(4)/time*length);
lamp3 = linspace(0,1,times(5)/time*length);
step3 = ones(1,times(6)/time*length);
lamp4 = linspace(1,0,times(7)/time*length);
step4 = zeros(1,times(8)/time*length);
u = [u1*[lamp1 step1 lamp2 step2] u2*[lamp3 step3 lamp4 step4]];


figure;
for i = 1:(length)/robo.const.ratiovel
    %plot trajectory
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    hold on;
    %plot robot posture
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.const.ratiovel
    robo = robo.shiftx;
    robo.u = u(:,robo.const.ratiovel*(i-1)+j);
    end
end
%%
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
figure;
for i = 1:3
subplot(3,1,i);
plot(robo.Tlog,robo.Xlog(3+i,:));
xlim([0 time]);
end
