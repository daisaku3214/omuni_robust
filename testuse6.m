close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter

robo = omunirobot(parameter,Deltat,Matrix,[p0;zeros(3+ell,1)],zeros(ell,1));
Q = diag([1000 1000 5000 100 100 500 0.001*ones(1,ell)]);
R = diag([0.01*ones(1,ell)]);
Aall = [zeros(3,3) eye(3,3) zeros(3,ell);
       zeros(3+ell,3) A];
Ball = [zeros(3,ell);B];
Kd = lqrd(Aall,Ball,Q,R,dt_simvel);

ref = [1 1 pi/4]';
time = 2;
length = time/robo.Deltat.simvel;
figure;
for i = 1:(length)/robo.const.ratiovel
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    hold on;
    %plot robot posture
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.const.ratiovel
    robo = robo.shiftx;
    Ttrans = robo.Ttrans;
    xe =  Ttrans*(robo.x(1:3)-ref);
    ve =  Ttrans*robo.x(4:6);
    robo.u = -Kd*([xe;ve;robo.x(7:end)]);
    
    end
end
%%
makemoviefromomunirobot_withref(robo,limit,ref,'movieK22deg');
%%
%plot grip force
hfin = figure;
set(hfin,'position',[150 150 1400 ell*150]);
Fmax = robo.para.mu_m'.*robo.const.N;
for i = 1:ell
subplot(ell,3,1+3*(i-1));
plot(robo.Tlog,robo.Flog(i,:),'-m','LineWidth',2);
xlim([0 time]);
hold on;grid on;
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[ Fmax(i)  Fmax(i)],'--g','LineWidth',2);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[-Fmax(i) -Fmax(i)],'--g','LineWidth',2);
ylabel(strcat(strcat('$$F_',num2str(i)),'[N]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,1);
title('Grip force', 'interpreter', 'latex');
%plot motor current
for i = 1:ell
subplot(ell,3,2+3*(i-1));
plot(robo.Tlog,robo.Xlog(6+i,:),'-m','LineWidth',2);
hold on;grid on;
xlim([0 time]);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[ robo.para.ilim_m(i)  robo.para.ilim_m(i)],'--g','LineWidth',2);
plot([robo.Tlog(1,1) robo.Tlog(1,end)],[-robo.para.ilim_m(i) -robo.para.ilim_m(i)],'--g','LineWidth',2);
ylabel(strcat(strcat('$$i_',num2str(i)),'[A]$$'), 'interpreter', 'latex');
end
xlabel('time [sec]', 'interpreter', 'latex');
subplot(ell,3,2);
title('Motor current', 'interpreter', 'latex');

%plot input voltage
for i = 1:ell
subplot(ell,3,3+3*(i-1));
plot(robo.Tlog,robo.Ulog(i,:),'-m','LineWidth',2);
hold on;grid on;
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
plot(robo.Tlog,robo.Xlog(i,:),'-m','LineWidth',2);
hold on; grid on;
plot([robo.Tlog(1,1) robo.Tlog(1,end)],ref(i)*[1 1],'--c','LineWidth',2)
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
plot(robo.Tlog,robo.Xlog(3+i,:),'-m','LineWidth',2);
hold on; grid on;
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
plot(robo.Tlog,robo.Xlog(i,:)-ref(i),'-m','LineWidth',2);
hold on; grid on;
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