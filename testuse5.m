close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter

parameter.rgast = -parameter.rgast;parameter.thetagast=pi/2+parameter.thetagast;parameter.I0 = 10;parameter.d_m =parameter.d_m*10; 
parameter.D_m = parameter.D_m*diag([1.01 1 0.99 1.02]);parameter.L_m = parameter.L_m*diag([0.9 1.1 1.2 0.8]);parameter.R_m = parameter.R_m*diag([2 1.2 1.4 0.8]);
parameter.m0 = 1.1*parameter.m0;parameter.J_m = parameter.J_m*diag([1.2 0.8 2 0.9]);parameter.eta_m = parameter.eta_m*diag([0.9 0.95 0.9 0.85]);

x0 = [p0+pe;zeros(3+ell,1)];
hatx0 = [p0;zeros(3+ell,1)];
robo = omunimachine(parameter,Deltat,Matrix,Uncertain,sensparas,controller_parameter,x0,hatx0,zeros(ell,1));
robo = robo.setKd(Kd);

mode = 'LQR_single_filt';

%input example:ramp input:move circle
track1_define;

figure;
for i = 1:(length)/robo.robot.const.ratiovel
    %plot trajectory
    plot(p(1,:),p(2,:),'--c');
    hold on;
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    %plot robot posture
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.robot.const.ratiovel
    xref = [p(:,robo.robot.const.ratiovel*(i-1)+j);
            v(:,robo.robot.const.ratiovel*(i-1)+j);
             zeros(ell,1)];
    u_ff = u(:,robo.robot.const.ratiovel*(i-1)+j);
    robo = robo.control_shift(xref,u_ff,mode);
    end
end
%%
figure;
makemoviefromomunirobot_withref(robo,limit,p,'movieK22deg');
%%
%plot input logger
robo.plotinputlogger(viewestimatei,viewnominalu);
%plot output logger
robo.plotoutputlogger(viewestimatex);
%plot estimate error
robo.plotestimateerrorlogger;