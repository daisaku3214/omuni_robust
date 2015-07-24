close all;
clear;
robot_define;           %definition robot parameter
control_define;         %definition contorol parameter

x0 = [p0;zeros(3+ell,1)];
hatx0 = [p0;zeros(3+ell,1)];
robo = omunimachine(parameter,Deltat,Matrix,Uncertain,sensparas,controller_parameter,x0,hatx0,zeros(ell,1));
robo = robo.setKd(Kdreg);

mode = 'LQR_single_filt';

ref = [1 1 pi/4]';
time = 2;
length = time/robo.robot.Deltat.simvel;
figure;
for i = 1:(length)/robo.robot.const.ratiovel
    plot(robo.Xlog(1,:),robo.Xlog(2,:),'m');
    hold on;
    %plot robot posture
    plot_omunirobot(robo,limit);
    drawnow;
    for j = 1:robo.robot.const.ratiovel
    xref = [ref;
            zeros(3,1);
            zeros(ell,1)];
    u_ff = zeros(ell,1);
    robo = robo.control_shift(xref,u_ff,mode);
    end
end
%%
figure;
makemoviefromomunirobot_withref(robo,limit,ref,'movieK22deg');
%%
%plot input logger
robo.plotinputlogger(viewestimatei,viewnominalu);
%plot output logger
robo.plotoutputlogger(viewestimatex);
%plot estimate error
robo.plotestimateerrorlogger;