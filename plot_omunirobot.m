function plot_omunirobot(omunirobot,limit)
tirewidth = omunirobot.robot.para.h_m;
dir = mean(omunirobot.robot.para.r)*[cos(omunirobot.robot.x(3));sin(omunirobot.robot.x(3))];
tiredir = [cos(omunirobot.robot.para.alpha + omunirobot.robot.x(3));...
           sin(omunirobot.robot.para.alpha + omunirobot.robot.x(3))];
Xbody = omunirobot.robot.para.r.*tiredir(1,:) + omunirobot.robot.x(1);
Ybody = omunirobot.robot.para.r.*tiredir(2,:) + omunirobot.robot.x(2);
fill(Xbody,Ybody,'w');
hold on;
Xtire = zeros(4,omunirobot.robot.para.ell);
Xtire(1,:) = (omunirobot.robot.para.r - tirewidth).*tiredir(1,:) + ...
              omunirobot.robot.para.D_m/2.*tiredir(2,:)+omunirobot.robot.x(1);
Xtire(2,:) = (omunirobot.robot.para.r + tirewidth).*tiredir(1,:) + ...
              omunirobot.robot.para.D_m/2.*tiredir(2,:)+omunirobot.robot.x(1);
Xtire(3,:) = (omunirobot.robot.para.r + tirewidth).*tiredir(1,:) - ...
              omunirobot.robot.para.D_m/2.*tiredir(2,:)+omunirobot.robot.x(1);
Xtire(4,:) = (omunirobot.robot.para.r - tirewidth).*tiredir(1,:) - ...
              omunirobot.robot.para.D_m/2.*tiredir(2,:)+omunirobot.robot.x(1);
Ytire = zeros(4,omunirobot.robot.para.ell);
Ytire(1,:) = (omunirobot.robot.para.r - tirewidth).*tiredir(2,:) - ...
              omunirobot.robot.para.D_m/2.*tiredir(1,:)+omunirobot.robot.x(2);
Ytire(2,:) = (omunirobot.robot.para.r + tirewidth).*tiredir(2,:) - ...
              omunirobot.robot.para.D_m/2.*tiredir(1,:)+omunirobot.robot.x(2);
Ytire(3,:) = (omunirobot.robot.para.r + tirewidth).*tiredir(2,:) + ...
              omunirobot.robot.para.D_m/2.*tiredir(1,:)+omunirobot.robot.x(2);
Ytire(4,:) = (omunirobot.robot.para.r - tirewidth).*tiredir(2,:) + ...
              omunirobot.robot.para.D_m/2.*tiredir(1,:)+omunirobot.robot.x(2);
for i = 1:omunirobot.robot.para.ell
 fill(Xtire(:,i),Ytire(:,i),'b');
 plot([omunirobot.robot.x(1);Xbody(:,i)],[omunirobot.robot.x(2);Ybody(:,i)],'k:');
end
plot([omunirobot.robot.x(1);2*dir(1,1)+omunirobot.robot.x(1)],...
     [omunirobot.robot.x(2);2*dir(2,1)+omunirobot.robot.x(2)],'m--');
plot([omunirobot.robot.x(1);2*omunirobot.robot.para.r(1)+omunirobot.robot.x(1)],...
     [omunirobot.robot.x(2);omunirobot.robot.x(2)],'k--');
plot(omunirobot.robot.x(1),omunirobot.robot.x(2),'mo','MarkerFaceColor','m');
robogast = [omunirobot.robot.x(1)+omunirobot.robot.para.rgast*...
            cos(omunirobot.robot.x(3)+omunirobot.robot.para.thetagast);...
            omunirobot.robot.x(2)+omunirobot.robot.para.rgast*...
            sin(omunirobot.robot.x(3)+omunirobot.robot.para.thetagast)];
plot(robogast(1),robogast(2),'go','MarkerFaceColor','g');
xlim(limit(1,:));
ylim(limit(2,:));
axis square;

hold off;
end