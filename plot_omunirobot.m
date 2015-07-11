function plot_omunirobot(omunirobot,limit)
tirewidth = omunirobot.para.h_m;
dir = mean(omunirobot.para.r)*[cos(omunirobot.x(3));sin(omunirobot.x(3))];
tiredir = [cos(omunirobot.para.alpha + omunirobot.x(3));...
           sin(omunirobot.para.alpha + omunirobot.x(3))];
Xbody = omunirobot.para.r.*tiredir(1,:) + omunirobot.x(1);
Ybody = omunirobot.para.r.*tiredir(2,:) + omunirobot.x(2);
fill(Xbody,Ybody,'w');
hold on;
Xtire = zeros(4,omunirobot.para.ell);
Xtire(1,:) = (omunirobot.para.r - tirewidth).*tiredir(1,:) + ...
              omunirobot.para.D_m/2.*tiredir(2,:)+omunirobot.x(1);
Xtire(2,:) = (omunirobot.para.r + tirewidth).*tiredir(1,:) + ...
              omunirobot.para.D_m/2.*tiredir(2,:)+omunirobot.x(1);
Xtire(3,:) = (omunirobot.para.r + tirewidth).*tiredir(1,:) - ...
              omunirobot.para.D_m/2.*tiredir(2,:)+omunirobot.x(1);
Xtire(4,:) = (omunirobot.para.r - tirewidth).*tiredir(1,:) - ...
              omunirobot.para.D_m/2.*tiredir(2,:)+omunirobot.x(1);
Ytire = zeros(4,omunirobot.para.ell);
Ytire(1,:) = (omunirobot.para.r - tirewidth).*tiredir(2,:) - ...
              omunirobot.para.D_m/2.*tiredir(1,:)+omunirobot.x(2);
Ytire(2,:) = (omunirobot.para.r + tirewidth).*tiredir(2,:) - ...
              omunirobot.para.D_m/2.*tiredir(1,:)+omunirobot.x(2);
Ytire(3,:) = (omunirobot.para.r + tirewidth).*tiredir(2,:) + ...
              omunirobot.para.D_m/2.*tiredir(1,:)+omunirobot.x(2);
Ytire(4,:) = (omunirobot.para.r - tirewidth).*tiredir(2,:) + ...
              omunirobot.para.D_m/2.*tiredir(1,:)+omunirobot.x(2);
for i = 1:omunirobot.para.ell
 fill(Xtire(:,i),Ytire(:,i),'b');
 plot([omunirobot.x(1);Xbody(:,i)],[omunirobot.x(2);Ybody(:,i)],'k:');
end
plot([omunirobot.x(1);2*dir(1,1)+omunirobot.x(1)],...
     [omunirobot.x(2);2*dir(2,1)+omunirobot.x(2)],'m--');
plot([omunirobot.x(1);2*omunirobot.para.r(1)+omunirobot.x(1)],...
     [omunirobot.x(2);omunirobot.x(2)],'k--');
plot(omunirobot.x(1),omunirobot.x(2),'mo','MarkerFaceColor','m');
robogast = [omunirobot.x(1)+omunirobot.para.rgast*...
            cos(omunirobot.x(3)+omunirobot.para.thetagast);...
            omunirobot.x(2)+omunirobot.para.rgast*...
            sin(omunirobot.x(3)+omunirobot.para.thetagast)];
plot(robogast(1),robogast(2),'go','MarkerFaceColor','g');
xlim(limit(1,:));
ylim(limit(2,:));
axis square;

hold off;
end