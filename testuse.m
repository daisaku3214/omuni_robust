robot_define;           %definition robot parameter
control_define;         %definition contorol parameter
time = 1;               %[sec] simuration time
length = time/Deltat.simvel;
robo = omunirobot(parameter,Deltat,[p0;zeros(3+ell,1)],zeros(ell,1));

u0 = [10 -10 -10 10]';
lamp1 = linspace(0,1,length/4);
step = ones(1,length/2);
lamp2 = linspace(1,0,length/4);
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