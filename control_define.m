w = logspace(-5,5,1000);%frequency of considering range
numP = 200;             %number of real plants
p0 = [0;0;0];           %initial posture of machine
limit = [-5 5; -5 5];   %limitation of field
Qv = diag([100000 100000 1000000 1*ones(1,ell)]);
Rv = diag([1*ones(1,ell)]);
Kvd = lqrd(A,B,Qv,Rv,dt_simvel);
Kv = lqr(A,B,Qv,Rv);
system = ss(A,B,eye(size(A)),zeros(size(B)));
dsystem = c2d(system,dt_simvel);
Lv = loopsens(system,Kv);
Lvd = loopsens(dsystem,Kvd);

Ap = [zeros(3,3) eye(3,3) zeros(3,ell);
       zeros(3+ell,3) Lv.To.a];
Bp = [zeros(3,3+ell);Lv.To.b];
Cp = eye(6+ell);
Dp = zeros(6+ell,3+ell);
psys = ss(Ap,Bp,Cp,Dp);
Qp2 = diag([10 10 20 0.00001 0.00001 0.000005 0.0000001*ones(1,ell)]);
Rp2 = diag([0.1 0.1 0.5 0.0000001*ones(1,ell)]);
Kpd = lqrd(Ap,Bp,Qp2,Rp2,dt_simvel);
Kp = lqr(Ap,Bp,Qp2,Rp2);
Lp = loopsens(psys,Kp);
Lpd= loopsens(psys,Kpd);

% Qp = diag([1000000 1000000 500000 1*ones(1,ell)]);
% Rp = diag([10*ones(1,3) 1*ones(1,ell)]);
% Kpd = lqrd(Lv.To.a,Lv.To.b,Qp,Rp,dt_simvel);
% Kp = lqr(Lv.To.a,Lv.To.b,Qp,Rp);
% Lp=loopsens(Lv.To,Kp);
% Lpcd = loopsens(c2d(Lv.To,dt_simvel),Kpd);
% Lpdd = loopsens(Lvd.To,Kpd);

