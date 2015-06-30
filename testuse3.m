close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
robo = omunirobot(parameter,Deltat,Matrix,[p0+diag([0.01 0.01 0.001])*randn(3,1);zeros(3+ell,1)],zeros(ell,1));
%input example:ramp input:move circle
% R = 10;
% omega = -pi/10;
% phi0 = -pi/6;
% x0 = R*cos(phi0);
% y0 = R*sin(phi0);
% theta0 = 0;
% angle = -4;
% tra = TRA_circle(R,x0,y0,theta0,pi+phi0,angle);
%[p,v,t] = calc_ref(tra,0,-pi/10,omega,omega,dt_simvel);

angle = pi/10;
phi1 = 0;
phi2 = 5*pi/4;
phi3 = pi/2;
se = 4;
vratio = 4;
accraito = 1;
se1 = se;
se2 = se*sqrt(2);
se3 = se;
vs1 = se1/vratio;
vs2 = se2/vratio;
vs3 = se3/vratio;

tra1 = TRA_line(0,0,0,phi1,angle);
tra2 = TRA_line(se1*cos(phi1),se1*sin(phi1),0,phi2,angle);
tra3 = TRA_line(se1*cos(phi1)+se2*cos(phi2),se1*sin(phi1)+se2*sin(phi2),0,phi3,angle);
[p1,v1,t1] = calc_ref(tra1,0,se1,vs1,vs1/accraito,dt_simvel);
[p2,v2,t2] = calc_ref(tra2,0,se2,vs2,vs2/accraito,dt_simvel);
[p3,v3,t3] = calc_ref(tra3,0,se3,vs3,vs3/accraito,dt_simvel);
p = [p1 p2 p3];
v = [v1 v2 v3];
t = [t1 t1(end)+t2 t1(end)+t2(end)+t3];

time = t(end);
length = size(t,2);
u = v2omega*[cos(p(3,1)) -sin(p(3,1)) 0;-sin(p(3,1)) cos(p(3,1)) 0; 0 0 1]*v(:,1);
for i = 1:length-1;
u = [u v2omega*[cos(p(3,i)) -sin(p(3,i)) 0;-sin(p(3,i)) cos(p(3,i)) 0; 0 0 1]*v(:,i)];
end
figure;
for i = 1:(length)/robo.const.ratiovel
    %plot trajectory
    plot(p(1,:),p(2,:),'--c');
    hold on;
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    %plot robot posture
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.const.ratiovel
    robo = robo.shiftx;
    xe = robo.x(1:3) - p(1:3,robo.const.ratiovel*(i-1)+j);
    vast = robo.Ttrans*v(:,robo.const.ratiovel*(i-1)+j);
    vref = [vast;zeros(ell,1)] + Kpd*([robo.Ttrans*xe;zeros(ell,1)]);

    robo.u = -Kvd*([robo.Ttrans*robo.x(4:6);robo.x(7:end)]-vref)...
             +v2omega*vast;

    Vrefworld(:,robo.const.ratiovel*(i-1)+j) =...
    [cos(p(3,robo.const.ratiovel*(i-1)+j)) -sin(p(3,robo.const.ratiovel*(i-1)+j)) 0;
    sin(p(3,robo.const.ratiovel*(i-1)+j)) cos(p(3,robo.const.ratiovel*(i-1)+j)) 0;
    0 0 1]*v(1:3,robo.const.ratiovel*(i-1)+j);
    end
end
%%
makemoviefromomunirobot_withref(robo,limit,p,'movieK22deg');
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
plot(t,u(i,:),'--m');
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
plot(t,p(i,:),'--m')
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
plot(t, v(i,:),'--m')
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
num1 = size(robo.Tlog,2);
num2 = size(p,2);
num = min(num1,num2);
%plot robot posture error
for i = 1:3
subplot(3,3,3+3*(i-1));
plot(robo.Tlog,robo.Xlog(i,1:num)-p(i,1:num),'-b');
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