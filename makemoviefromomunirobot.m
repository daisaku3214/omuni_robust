function makemoviefromomunirobot(omunirobot,limit,title)
mvlgth = size(omunirobot.Xlog,2);
writerObj = VideoWriter(title);
writerObj.FrameRate = 1/omunirobot.robot.Deltat.simpos;
open(writerObj);
for i = 1:mvlgth
  if (mod(i-1,omunirobot.const.ratiovel)==0)
  plot(omunirobot.Xlog(1,1:i),omunirobot.Xlog(2,1:i),'m','linewidth',2);
  hold on;
  omunirobot.robot.x = omunirobot.Xlog(:,i);
  plot_omunirobot(omunirobot,limit);
  drawnow;
  hold off;
  frame = getframe;
  writeVideo(writerObj,frame);  
  end
end

close(writerObj);
end