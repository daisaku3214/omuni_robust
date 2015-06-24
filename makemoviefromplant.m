function makemoviefromplant(plants,limit,title)
mvlgth = size(plants.nominal.Xlog,2);
writerObj = VideoWriter(title);
writerObj.FrameRate = 1/plants.nominal.Deltat.simpos;
open(writerObj);
for i = 1:mvlgth
  if (mod(i-1,plants.nominal.const.ratiovel)==0)
  plot(plants.nominal.Xlog(1,1:i),plants.nominal.Xlog(2,1:i),'m','linewidth',2);
  hold on;
  for j = 1:plants.numP
    plot(plants.real(j).Xlog(1,1:i),plants.real(j).Xlog(2,1:i),'b');
  end
  plot(plants.nominal.Xlog(1,1:i),plants.nominal.Xlog(2,1:i),'m','linewidth',2);
  for j = 1:plants.numP
    hold on;
    plants.real(j).x = plants.real(j).Xlog(:,i);
    plot_omunirobotmini(plants.real(j),limit);
    hold on;
  end
  plants.nominal.x = plants.nominal.Xlog(:,i);
  plot_omunirobot(plants.nominal,limit);
  drawnow;
  hold off;
  frame = getframe;
  writeVideo(writerObj,frame);
  end
end
end