function plot_plants(omuniplants,limit)
plot(omuniplants.nominal.Xlog(1,:),omuniplants.nominal.Xlog(2,:),'m');
hold on;
for i = 1:omuniplants.numP
    plot(omuniplants.real(i).Xlog(1,:),omuniplants.real(i).Xlog(2,:),'b');
end
plot(omuniplants.nominal.Xlog(1,:),omuniplants.nominal.Xlog(2,:),'m');
plot_omunirobot(omuniplants.nominal,limit);
end