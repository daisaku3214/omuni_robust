close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
robo = omunirobot(parameter,Deltat,Matrix,[p0+diag([0.01 0.01 0.001])*randn(3,1);zeros(3+ell,1)],zeros(ell,1));
%input example:ramp input:move circle
 R = 2;
 omega = -2*pi/8;
 phi0 = 2*pi/3;
 x0 = R*cos(phi0);
 y0 = R*sin(phi0);
 se0 = -2*pi;
 theta0 = 0;
 angle0 = -0.5;
 tra = TRA_circle(R,x0,y0,theta0,pi+phi0,angle0);
[p0,v0,t0] = calc_ref(tra,0,se0,omega,omega,dt_simvel);

angle1 = pi/8;
angle2 = 0;
angle3 = -pi/8;
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

tra1 = TRA_line(p0(1,end),p0(2,end),p0(3,end),phi1,angle1);
[p1,v1,t1] = calc_ref(tra1,0,se1,vs1,vs1/accraito,dt_simvel);
t1 = t1 + t0(end)+dt_simvel;
tra2 = TRA_line(p1(1,end),p1(2,end),p1(3,end),phi2,angle2);
[p2,v2,t2] = calc_ref(tra2,0,se2,vs2,vs2/accraito,dt_simvel);
t2 = t2 + t1(end)+dt_simvel;
tra3 = TRA_line(p2(1,end),p2(2,end),p2(3,end),phi3,angle3);
[p3,v3,t3] = calc_ref(tra3,0,se3,vs3,vs3/accraito,dt_simvel);
t3 = t3 + t2(end)+dt_simvel;

p = [p0 p1 p2 p3];
v = [v0 v1 v2 v3];
t = [t0 t1 t2 t3];

time = t(end);
length = size(t,2);
vm = zeros(size(v));
for i = 1:length-1;
vm(:,i) = [cos(p(3,i)) sin(p(3,i)) 0;-sin(p(3,i)) cos(p(3,i)) 0; 0 0 1]*v(:,i);
end
u = v2omega*vm;
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
    vref = [vast;zeros(ell,1)] + Kpd*([robo.Ttrans*xe;robo.x(7:end,1)]);

    robo.u = -Kvd*([robo.Ttrans*robo.x(4:6);robo.x(7:end)]-vref)...
             +v2omega*vast;

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