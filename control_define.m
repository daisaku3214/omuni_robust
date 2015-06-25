w = logspace(-3,5,1000);%frequency of considering range
numP = 200;             %number of real plants
p0 = [0;0;0];           %initial posture of machine
limit = [-5 5; -5 5];   %limitation of field
Q = diag([10000 10000 10000 1*ones(1,ell)]);
R = diag([1*ones(1,ell)]);
Kd = lqrd(A,B,Q,R,dt_simvel);
K = lqr(A,B,Q,R);