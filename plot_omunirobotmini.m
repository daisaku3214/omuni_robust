function plot_omunirobotmini(omunirobot,limit)
dir = mean(omunirobot.para.r)*[cos(omunirobot.x(3));sin(omunirobot.x(3))];

plot([omunirobot.x(1);2*dir(1,1)+omunirobot.x(1)],...
     [omunirobot.x(2);2*dir(2,1)+omunirobot.x(2)],'b--');
hold on;
plot([omunirobot.x(1);2*omunirobot.para.r(1)+omunirobot.x(1)],...
     [omunirobot.x(2);omunirobot.x(2)],'k--');
plot(omunirobot.x(1),omunirobot.x(2),'bo','MarkerFaceColor','b');
xlim(limit(1,:));
ylim(limit(2,:));
axis square;

hold off;
end