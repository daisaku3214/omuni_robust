w = logspace(-5,5,1000);%frequency of considering range
numP = 200;             %number of real plants
p0 = [0;0;0];           %initial posture of machine
pe = [0.5;0.0;0.0];   %initial posture error
limit = [-5 5; -5 5];   %limitation of field
viewestimatex = 1;      %plot estimate state value in log plot
viewestimatei = 0;      %plot estimate state value in log plot
viewnominalu = 0;       %plot nominal input value without saturation
%%
%define of state estimator
alfp = 0.05;
alfENCvsJRO = 0.0;
fpxy = 2*pi*15;
fvxy = 2*pi*15;
fro= 2*pi*10;
fu = 2*pi*25;
%good fpxy  = 15[Hz]
%good fvxy  = 15[Hz]
%good fro   =  5[Hz]
%good fu    = 25[Hz]
%higher fpxy = 100[Hz]
%higher fvxy = 100[Hz]
%higher fro  = 15[Hz]
%higher fu   = 60[Hz]
%lower fpxy = 2.5[Hz]
%lower fvxy = 2.5[Hz]
%lower fro  = 2.5[Hz]
%lower fu    = 15[Hz]

controller_parameter = struct('alfp',alfp,'alfENCvsJRO',alfENCvsJRO,'fpxy',fpxy,...
                       'fvxy',fvxy,'fro',fro,'fu',fu,'LRFratio',LRF_dt/dt_simvel);

%%
% Qv = diag([100000 100000 1000000 0.01*ones(1,ell)]);
% Rv = diag([1*ones(1,ell)]);
% Kvd = lqrd(A,B,Qv,Rv,dt_simvel);
% Kv = lqr(A,B,Qv,Rv);
% Kvdreg = [Kvd(:,1:3) zeros(ell,ell)];
% Kvreg = [Kv(:,1:3) zeros(ell,ell)];
% system = ss(A,B,eye(size(A)),zeros(size(B)));
% dsystem = c2d(system,dt_simvel);
% Lv = loopsens(system,Kv);
% Lvd = loopsens(dsystem,Kvd);
% Lvreg = loopsens(system,Kvreg);
% Lvdreg = loopsens(dsystem,Kvdreg);
% 
% Ap = [zeros(3,3) eye(3,3) zeros(3,ell);
%        zeros(3+ell,3) Lv.To.a];
% Bp = [zeros(3,3+ell);Lv.To.b];
% Cp = eye(6+ell);
% Dp = zeros(6+ell,3+ell);
% psys = ss(Ap,Bp,Cp,Dp);
% Qp = diag([10 10 20 0.0001 0.0001 0.00005 0.0000001*ones(1,ell)]);
% Rp = diag([1 1 5 0.0000001*ones(1,ell)]);
% Kpd = lqrd(Ap,Bp,Qp,Rp,dt_simvel);
% Kp = lqr(Ap,Bp,Qp,Rp);
% Kpdreg=[Kpd(1:3,1:6) zeros(3,ell);zeros(ell,6+ell)];
% Kpreg=[Kp(1:3,1:6) zeros(3,ell);zeros(ell,6+ell)];
% Lp = loopsens(psys,Kp);
% Lpd= loopsens(psys,Kpd);
% Lpreg = loopsens(psys,Kpreg);
% Lpdreg= loopsens(psys,Kpdreg);

Q = diag([10000 10000 50000 1000 1000 5000 0.001*ones(1,ell)]);
R = diag([0.01*ones(1,ell)]);
Aall = [zeros(3,3) eye(3,3) zeros(3,ell);
       zeros(3+ell,3) A];
Ball = [zeros(3,ell);B];
Call = eye(6+ell,6+ell);
Dall = zeros(6+ell,ell);
allsys = ss(Aall,Ball,Call,Dall);
alldsys = c2d(allsys,Deltat.simvel);

Kd = lqrd(Aall,Ball,Q,R,dt_simvel);
Kdreg = [Kd(:,1:6) zeros(ell,ell)];
Ldall = loopsens(alldsys,Kd);
Ldallreg = loopsens(alldsys,Kdreg);


% Zirrep = care(A',C',B*B');
% Hirrep = -Zirrep*C';
% Nirrep = ss(A+Hirrep*C,B,C,zeros(size(C,1),size(B,2)));
% Mirrep = ss(A+Hirrep*C,Hirrep,C,eye(size(C,1),size(Hirrep,2)));

