close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
robo = omunirobot(parameter,Deltat,Matrix,[p0+[0;0.01;0.01];zeros(3+ell,1)],zeros(ell,1));
tsta = 1.6820;   %[sec]
stanum = tsta/Deltat.simvel;
Q = diag([1000 1000 5000 100 100 500 0.001*ones(1,ell)]);
R = diag([0.01*ones(1,ell)]);
Aall = [zeros(3,3) eye(3,3) zeros(3,ell);
       zeros(3+ell,3) A];
Ball = [zeros(3,ell);B];
Kd = lqrd(Aall,Ball,Q,R,dt_simvel);

%input example:ramp input:move circle
R = 2;
omega = -2*pi/7.5;
acc_ratio = 2;
phi0 = 2*pi/3;
x0 = R*cos(phi0);
y0 = R*sin(phi0);
se0 = -2*pi;
theta0 = 0;
angle0 = -0.25;
tra = TRA_circle(R,x0,y0,theta0,pi+phi0,angle0);
[p0,v0,t0] = calc_ref(tra,0,se0,omega,omega/acc_ratio,dt_simvel);

angle1 = -pi/16;
angle2 = 0;
angle3 = -pi/16;
phi1 = 0;
phi2 = 5*pi/4;
phi3 = pi/2;
se = 4;
vratio = 2;
accraito = 1;
se1 = se;
se2 = se*sqrt(2);
se3 = se;
vs1 = se1/vratio;
vs2 = 1.1*se2/vratio;
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

p4 = zeros(3,1)*ones(1,stanum);
v4 = v3(:,end)*ones(1,stanum);
t4 = t3(:,end)+dt_simvel+linspace(0,tsta,stanum);

p = [p0 p1 p2 p3 p4];
v = [v0 v1 v2 v3 v4];
t = [t0 t1 t2 t3 t4];

time = t(end);
length = size(t,2);
vm = zeros(size(v));
u = omega2Vol*v2omega*vm;
for i = 1:length-1;
vm(:,i) = [cos(p(3,i)) sin(p(3,i)) 0;-sin(p(3,i)) cos(p(3,i)) 0; 0 0 1]*v(:,i);
u(:,i) = pinv(Matrix.dB)*([vm(:,i+1);zeros(ell,1)]-Matrix.dA*[vm(:,i);zeros(ell,1)]);
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
    Ttrans = robo.Ttrans;
    vast = Ttrans*v(:,robo.const.ratiovel*(i-1)+j);
    xe =  Ttrans*(robo.x(1:3) - p(1:3,robo.const.ratiovel*(i-1)+j));
    ve =  Ttrans*robo.x(4:6) - vast;
    robo.u = -Kd*([xe;ve;robo.x(7:end)])+u(:,robo.const.ratiovel*(i-1)+j);
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
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'k')  %orgin line
hold on; grid on;
plot(robo.Tlog,robo.Flog(i,:),'-m','LineWidth',2);
xlim([0 time]);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[ Fmax(i)  Fmax(i)],':g','LineWidth',2);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[-Fmax(i) -Fmax(i)],':g','LineWidth',2);
ylabel(strcat(strcat('$$F_',num2str(i)),'[N]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,1);
title('Grip force', 'interpreter', 'latex');
%plot motor current
for i = 1:ell
subplot(ell,3,2+3*(i-1));
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'k')  %orgin line
hold on; grid on;
plot(robo.Tlog,robo.Xlog(6+i,:),'-m','LineWidth',2);
xlim([0 time]);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[ robo.para.ilim_m(i)  robo.para.ilim_m(i)],':g','LineWidth',2);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[-robo.para.ilim_m(i) -robo.para.ilim_m(i)],':g','LineWidth',2);
ylabel(strcat(strcat('$$i_',num2str(i)),'[A]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,2);
title('Motor current', 'interpreter', 'latex');

%plot input voltage
for i = 1:ell
subplot(ell,3,3+3*(i-1));
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'k')  %orgin line
hold on; grid on;
plot(robo.Tlog,robo.Ulog(i,:),'-m','LineWidth',2);
plot(t,u(i,:),'--c','LineWidth',2);
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
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'k')  %orgin line
hold on; grid on;
plot(robo.Tlog,robo.Xlog(i,:),'-m','LineWidth',2);
plot(t,p(i,:),'--c','LineWidth',2)
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
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'k')  %orgin line
hold on; grid on;
plot(robo.Tlog,robo.Xlog(3+i,:),'-m','LineWidth',2);
plot(t, v(i,:),'--c','LineWidth',2)
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
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'k')  %orgin line
hold on; grid on;
plot(robo.Tlog,robo.Xlog(i,1:num)-p(i,1:num),'-m','LineWidth',2);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],zeros(1,2),'--c','LineWidth',2)
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