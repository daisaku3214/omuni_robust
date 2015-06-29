close all;
clear all;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
robo = omunirobot(parameter,Deltat,Matrix,[p0+diag([0.01 0.01 0.001])*randn(3,1);zeros(3+ell,1)],zeros(ell,1));
%input example:ramp input:move circle
R = 1;
omega = 2*pi*1/8;
phi0 = 1*pi/6;
x0 = -R*cos(phi0);
y0 = -R*sin(phi0);
theta0 = 0;

times = [1 7 1 1];     %[sec] the time of input shaping
time = sum(times);      %[sec] simuration time
length = time/Deltat.simvel;
lamp1 = linspace(0,1,times(1)/time*length);
step1 = ones(1,times(2)/time*length);
lamp2 = linspace(1,0,times(3)/time*length);
step2 = zeros(1,times(4)/time*length);
tlamp1 = linspace(0,times(1),times(1)/time*length);
theta1 = 0.5*omega/(times(1)).*(tlamp1-0).^2 + theta0;
x1 = x0 + R*cos(theta1+phi0);
y1 = y0 + R*sin(theta1+phi0);
tstep1 = linspace(times(1)+Deltat.simvel,times(1)+times(2),times(2)/time*length);
theta2 = omega*(tstep1-times(1))+ theta1(end);
x2 = x0 + R*cos(theta2+phi0);
y2 = y0 + R*sin(theta2+phi0);
tlamp2 = linspace(times(1)+times(2)+Deltat.simvel,times(1)+times(2)+times(3),times(3)/time*length);
theta3 = omega.*(tlamp2-times(2)-times(1))-0.5*omega/(times(3)).*(tlamp2-times(2)-times(1)).^2 + theta2(end);
x3 = x0 + R*cos(theta3+phi0);
y3 = y0 + R*sin(theta3+phi0);
tstep2 = linspace(times(1)+times(2)+times(3)+Deltat.simvel,times(1)+times(2)+times(3)+times(4),times(4)/time*length);
theta4 = theta3(end)*ones(1,times(4)/time*length);
x4 = x0 + R*cos(theta4+phi0);
y4 = y0 + R*sin(theta4+phi0);

x = [x1 x2 x3 x4;
     y1 y2 y3 y4;
     theta1 theta2 theta3 theta4;
     zeros(ell,length)];
v1 = [-omega*R*sin(phi0) omega*R*cos(phi0) omega]';
u1 = v2omega*v1; u = u1*[lamp1 step1 lamp2 step2];
u = [u u(:,end)];
v = v1*[lamp1 step1 lamp2 step2];
v = [v;zeros(ell,length)];
Vrefworld = zeros(3,length);
x = [x x(:,end)];

% x = [x1 x2 x3 x4;
%      y1 y2 y3 y4;
%      zeros(1,length);
%      zeros(ell,length)];
% v1 = [-omega*R*sin(phi0) omega*R*cos(phi0) 0]';
% v = v1*[lamp1 step1 lamp2 step2];
% v = robo.Tarray([theta1 theta2 theta3 theta4],v);
% u = v2omega*v;
% u = [u u(:,end)];
% v = [v;zeros(ell,length)];
% Vrefworld = zeros(3,length);
% x = [x x(:,end)];

%input example:ramp input:move linear
% times = [1 1 1 1 1 1 1 1];     %[sec] the time of input shaping
% time = sum(times);      %[sec] simuration time
% length = time/Deltat.simvel;
% u1 = v2omega*[2 0 0]';
% u2 = v2omega*[-2 -2 0]';
% v1 = [2 0 0 zeros(1,ell)]';
% v2 = [-2 -2 0 zeros(1,ell)]';
% lamp1 = linspace(0,1,times(1)/time*length);
% step1 = ones(1,times(2)/time*length);
% lamp2 = linspace(1,0,times(3)/time*length);
% step2 = zeros(1,times(4)/time*length);
% lamp3 = linspace(0,1,times(5)/time*length);
% step3 = ones(1,times(6)/time*length);
% lamp4 = linspace(1,0,times(7)/time*length);
% step4 = zeros(1,times(8)/time*length);
% u = [u1*[lamp1 step1 lamp2 step2] u2*[lamp3 step3 lamp4 step4]];
% v = [v1*[lamp1 step1 lamp2 step2] v2*[lamp3 step3 lamp4 step4]];

figure;
for i = 1:(length)/robo.const.ratiovel
    %plot trajectory
    plot(x(1,:),x(2,:),'--c');
    hold on;
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    %plot robot posture
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.const.ratiovel
    robo = robo.shiftx;
    xe = robo.x(1:3) - x(1:3,robo.const.ratiovel*(i-1)+j);
    vref = v(:,robo.const.ratiovel*(i-1)+j) + Kpd*([robo.Ttrans*xe;zeros(ell,1)]);

    robo.u = -Kvd*([robo.Ttrans*robo.x(4:6);robo.x(7:end)]-vref)...
             +u(:,robo.const.ratiovel*(i-1)+j);
    Vrefworld(:,robo.const.ratiovel*(i-1)+j) =...
    [cos(x(3,robo.const.ratiovel*(i-1)+j)) -sin(x(3,robo.const.ratiovel*(i-1)+j)) 0;
    sin(x(3,robo.const.ratiovel*(i-1)+j)) cos(x(3,robo.const.ratiovel*(i-1)+j)) 0;
    0 0 1]*v(1:3,robo.const.ratiovel*(i-1)+j);
    end
end
%%
makemoviefromomunirobot_withref(robo,limit,x,'movieK22deg');
%%
%plot grip force
hfin = figure;
set(hfin,'position',[150 150 1400 ell*150]);
Fmax = robo.para.mu_m'.*robo.const.N;
for i = 1:ell
subplot(ell,3,1+3*(i-1));
plot(robo.Tlog,robo.Flog(i,:),'-b');
xlim([0 time]);
hold on;grid on;
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[ Fmax(i)  Fmax(i)],'--c');
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[-Fmax(i) -Fmax(i)],'--c');
ylabel(strcat(strcat('$$F_',num2str(i)),'[N]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,1);
title('Grip force', 'interpreter', 'latex');
%plot motor current
for i = 1:ell
subplot(ell,3,2+3*(i-1));
plot(robo.Tlog,robo.Xlog(6+i,:),'-b');
hold on;grid on;
xlim([0 time]);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[ robo.para.ilim_m(i)  robo.para.ilim_m(i)],'--c');
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[-robo.para.ilim_m(i) -robo.para.ilim_m(i)],'--c');
ylabel(strcat(strcat('$$i_',num2str(i)),'[A]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,2);
title('Motor current', 'interpreter', 'latex');

%plot input voltage
for i = 1:ell
subplot(ell,3,3+3*(i-1));
plot(robo.Tlog,robo.Ulog(i,:),'-b');
hold on;grid on;
plot(robo.Tlog,u(i,:),'--m');
xlim([0 time]);
ylabel(strcat(strcat('$$V_',num2str(i)),'[V]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,3);
title('Input voltage', 'interpreter', 'latex');

%plot robot posture
hfout = figure;
set(hfout,'position',[150 150 1400 600]);
for i = 1:3
subplot(3,3,1+3*(i-1));
plot(robo.Tlog,robo.Xlog(i,:),'-b');
hold on; grid on;
plot(robo.Tlog,x(i,:),'--m')
xlim([0 time]);
    if i==1
        title('Robot posture', 'interpreter', 'latex');
        ylabel('$$x {[m]}$$', 'interpreter', 'latex');
        ylim(limit(1,:));
    end
    if i==2
        ylabel('$$y {[m]}$$', 'interpreter', 'latex');
        ylim(limit(2,:));
    end
    if i==3
        ylabel('$$\theta {[rad]}$$', 'interpreter', 'latex');
    end
end
xlabel('time [sec]');

%plot robot velosity at world
for i = 1:3
subplot(3,3,2+3*(i-1));
plot(robo.Tlog,robo.Xlog(3+i,:),'-b');
hold on; grid on;
plot(robo.Tlog,[0 Vrefworld(i,:)],'--m')
xlim([0 time]);
    if i==1
        title('Robot velosity at world', 'interpreter', 'latex');
        ylabel('$$\dot x {[m/s]}$$', 'interpreter', 'latex');
    end
    if i==2
        ylabel('$$\dot y {[m/s]}$$', 'interpreter', 'latex');
    end
    if i==3
        ylabel('$$\dot \theta {[rad/s]}$$', 'interpreter', 'latex');
    end
end
xlabel('time [sec]');

%plot robot posture error
for i = 1:3
subplot(3,3,3+3*(i-1));
plot(robo.Tlog,robo.Xlog(i,:)-x(i,:),'-b');
hold on; grid on;
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'--m')
xlim([0 time]);
    if i==1
        title('Robot posture error', 'interpreter', 'latex');
        ylabel('$$e_x {[m]}$$', 'interpreter', 'latex');
    end
    if i==2
        ylabel('$$e_y {[m]}$$', 'interpreter', 'latex');
    end
    if i==3
        ylabel('$$e_\theta {[rad]}$$', 'interpreter', 'latex');
    end
end
xlabel('time [sec]');