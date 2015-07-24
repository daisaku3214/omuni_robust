%input example:ramp input:move circle and linear
tsta = 1.6800;   %[sec]
stanum = tsta/Deltat.simvel;

R = 2;
omega = -2*pi/7.5;
acc_ratio = 2;
phi0 = 2*pi/3;
x0 = R*cos(phi0);
y0 = R*sin(phi0);
se0 = -2*pi;
theta0 = 0;
angle0 = -0.25;
tra = TRA_circle(R,x0,y0,theta0,pi+phi0,angle0);
[p0,v0,t0] = calc_ref(tra,0,se0,omega,omega/acc_ratio,dt_simvel);

angle1 = -pi/16;
angle2 = 0;
angle3 = -pi/16;
phi1 = 0;
phi2 = 5*pi/4;
phi3 = pi/2;
se = 4;
vratio = 2;
accraito = 2;
se1 = se;
se2 = se*sqrt(2);
se3 = se;
vs1 = se1/vratio;
vs2 = 1*se2/vratio;
vs3 = se3/vratio;

tra1 = TRA_line(p0(1,end),p0(2,end),p0(3,end),phi1,angle1);
[p1,v1,t1] = calc_ref(tra1,0,se1,vs1,vs1/accraito,dt_simvel);
t1 = t1 + t0(end)+dt_simvel;
tra2 = TRA_line(p1(1,end),p1(2,end),p1(3,end),phi2,angle2);
[p2,v2,t2] = calc_ref(tra2,0,se2,vs2,vs2/accraito,dt_simvel);
t2 = t2 + t1(end)+dt_simvel;
tra3 = TRA_line(p2(1,end),p2(2,end),p2(3,end),phi3,angle3);
[p3,v3,t3] = calc_ref(tra3,0,se3,vs3,vs3/accraito,dt_simvel);
t3 = t3 + t2(end)+dt_simvel;

p4 = zeros(3,1)*ones(1,stanum);
v4 = v3(:,end)*ones(1,stanum);
t4 = t3(:,end)+dt_simvel+linspace(0,tsta,stanum);

p = [p0 p1 p2 p3 p4];
v = [v0 v1 v2 v3 v4];
t = [t0 t1 t2 t3 t4];

time = t(end);
length = size(t,2);

vm = zeros(size(v));
u = omega2Vol*v2omega*vm;