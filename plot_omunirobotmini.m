function plot_omunirobotmini(omunirobot,limit)
dir = mean(omunirobot.robot.para.r)*[cos(omunirobot.robot.x(3));sin(omunirobot.robot.x(3))];

plot([omunirobot.robot.x(1);2*dir(1,1)+omunirobot.robot.x(1)],...
     [omunirobot.robot.x(2);2*dir(2,1)+omunirobot.robot.x(2)],'b--');
hold on;
plot([omunirobot.robot.x(1);2*omunirobot.robot.para.r(1)+omunirobot.robot.x(1)],...
     [omunirobot.robot.x(2);omunirobot.robot.x(2)],'k--');
plot(omunirobot.robot.x(1),omunirobot.robot.x(2),'bo','MarkerFaceColor','b');
xlim(limit(1,:));
ylim(limit(2,:));
axis square;

hold off;
end