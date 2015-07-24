function hatvw = v_rev(vw,dt)
hatvw = zeros(size(vw));
det = 1+(vw(3,1)*dt^2/2)^2;
err = vw(3,1)*dt^2/2;

hatvw(1,1) = vw(1,1)/det -err*vw(2,1)/det;
hatvw(2,1) = err*vw(1,1)/det + vw(2,1)/det;
hatvw(3,1) = vw(3,1);
end