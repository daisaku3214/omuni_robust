%%
%robot_define
%this file is definition of robot parameter
%%
%common definition
dt_simcir =1e-4;%sampling time for circuit equation
dt_simvel =1e-3;%sampling time for velocity simulation
dt_simpos =1e-2;%sampling time for position simulation
g = 9.80665;    %[m/s^2] the Acceleration of gravity
ell = 4;        %the number of motor
I0 = 8;        %[kgm^2]    z-axis inertia of robot other than motor units
m0 = 4;        %[kg]       mass of robot other than motor units
rgast = 0;      %[m]        the radius between COM of robot and robot origin
thetagast = 0;  %[rad]      the radian from x ast axis
alpha = linspace(0,2*pi-2*pi/ell,ell);
rc = 0.5;       %[m]        the radius between omni wheel and robot origin
r = rc*ones(1,ell);
rprim = sqrt(rgast^2.+r.^2 -2*rgast.*r.*cos(thetagast-alpha));
cbeta = (rprim./r + r./rprim -rgast.^2./(r.*rprim))/2;
numi = 0;       %the number of current sensor

%%
%motor unit definition
%at this mode, all motor is same
%this patteren, I will use RE35
Rc = 0.583;     %[ohm]      The internal resistance of the nominal motor
Lc = 0.191/1000;%[H]        The internal reactance of the nominal motor
KTc = 29.2/1000;%[Nm/A]     The Torque constant of the nominal motor
ROT = 328;      %[rpm/V]    The Rotational constant of the nominal motor
ROT=2*pi*ROT/60;%[rad/Vs]   The Rotational constant of the nominal motor
KEc = 1/(ROT);  %[Vs/rad]   The Back electromotive force constanto of the nominal motor
JAc = 79.2;     %[gcm^2]    The inertia of the rotor of the nominal motor
JAc =JAc/10^(7);%[kgm^2]    The inertia of the rotor of the nominal motor
Gc =  51;       %           The ratio of gear
etac = 0.7;     %           The efficiency of gear
muc = 0.8;      %           The Coefficient of static friction
ilimc = 6;     %[A]        the limit of motor current
Vlimc = 22;     %[V]        the limit of input voltage
%this pattern, I use omni wheel of Diameter 200mm
Dc = 0.4;       %[m]        The dirmeter of omni wheel
Jc = 1;      %[kgm^2]    The inertia of the omni whell
mc = 1.5;       %[kg]       The mass of a motor unit
dc = 0.01;       %[Nm/sec]    The declease ratio

R_m = Rc*ones(1,ell);
L_m = Lc*ones(1,ell);
KT_m = KTc*ones(1,ell);
KE_m = KEc*ones(1,ell);
JA_m = JAc*ones(1,ell);
G_m = Gc*ones(1,ell);
eta_m = etac*ones(1,ell);
D_m = Dc*ones(1,ell);
J_m = Jc*ones(1,ell);
m_m = mc*ones(1,ell);
d_m = dc*ones(1,ell);
mu_m = muc*ones(1,ell);
ilim_m=ilimc*ones(1,ell);
Vlim_m=Vlimc*ones(1,ell);

M = m0 + sum(m_m);
I = I0 + sum(rprim.^2.*m_m);
%%
sinalpha = sin(alpha);
cosalpha = cos(alpha);
a_m = 2*eta_m.*G_m.*KT_m./D_m;
b_m = 2*(J_m+eta_m.*G_m.^2.*JA_m)./D_m;
c_m = 2*b_m./D_m;
dD_m =2*d_m./D_m;
mMCss = sum(c_m.*sinalpha.^2)/M;
mMCcc = sum(c_m.*cosalpha.^2)/M;
mMCsc = sum(c_m.*sinalpha.*cosalpha)/M;
mMCrs = sum(r.*c_m.*sinalpha)/M;
mMCrc = sum(r.*c_m.*cosalpha)/M;
mMDss = sum(dD_m.*sinalpha.^2)/M;
mMDcc = sum(dD_m.*cosalpha.^2)/M;
mMDsc = sum(dD_m.*sinalpha.*cosalpha)/M;
mMDrs = sum(r.*dD_m.*sinalpha)/M;
mMDrc = sum(r.*dD_m.*cosalpha)/M;
mICbs = sum(rprim.*c_m.*cbeta.*sinalpha)/I;
mICbc = sum(rprim.*c_m.*cbeta.*cosalpha)/I;
mICbr = sum(rprim.*r.*c_m.*cbeta)/I;
mIDbs = sum(rprim.*dD_m.*cbeta.*sinalpha)/I;
mIDbc = sum(rprim.*dD_m.*cbeta.*cosalpha)/I;
mIDbr = sum(rprim.*dD_m.*cbeta)/I;
mMAs = a_m.*sinalpha/M;
mMAc = a_m.*cosalpha/M;
mIAr = rprim.*a_m.*cbeta/I;
mM = [1-mMCss mMCsc mMCrs-rgast*sin(thetagast);
      mMCsc 1-mMCcc -mMCrc+rgast*cos(thetagast);
      mICbs -mICbc 1-mICbr];
mD = [ mMDss -mMDsc -mMDrs;
     -mMDsc  mMDcc  mMDrc;
     -mIDbs  mIDbc  mIDbr];
mA = [-mMAs;mMAc;mIAr];
mAA = mM\mA;
mAD = mM\mD;
mR  = diag(R_m./L_m);
%wR  = mR;
mB2 = diag(ones(1,ell)./L_m);
%wB2 = mB2;
Komega = 2*G_m.*KE_m./(L_m.*D_m);
mKx = Komega.*sin(alpha);
mKy = Komega.*cos(alpha);
mKq = Komega.*r;
mK = [-mKx' mKy' mKq'];

A = [mAD mAA;
      mK -mR];
B = [zeros(3,ell);mB2];

C1 = [0 0 1 zeros(1,ell);
       (2*sinalpha./D_m)' -(2*cosalpha./D_m)' -(2*rc./D_m)' zeros(ell,ell)];
C2 = [zeros(numi,3) eye(numi,ell)];
C = [C1;C2];
D = zeros(size(C,1),ell);
v2omega = C1(2:ell+1,1:3);
omega2v = pinv(v2omega);


%Determining whether controllable
disp(size(A));
[~,~,~,~,k] = ctrbf(A,B,C);
disp(sum(k));
%Determining whether observability
[~,~,~,~,k] = obsvf(A,B,C);
disp(sum(k));
%Determine whether the detectability
lambda = eig(A);
[~,indextemp] = sort(real(lambda),'descend');
lambda = lambda(indextemp);
clear indextemp;
lambdatemp = 0;
for i = 1:size(A,1);
if (lambda(i) > 0)
    lambdatemp = lambdatemp + lambda(i);
end
end
lambda = lambdatemp;
clear lambdatemp;
disp(rank([C;A-lambda]));
clear lambda i k;


parameter = struct('g',g,'ell',ell,'I0',I0,'m0',m0,'rgast',rgast,...
                   'thetagast',thetagast,'alpha',alpha,'r',r,...
                   'numi',numi,'R_m',R_m,'L_m',L_m,'KT_m',KT_m,...
                   'KE_m',KE_m,'JA_m',JA_m,'G_m',G_m,'eta_m',eta_m,...
                   'D_m',D_m,'J_m',J_m,'m_m',m_m,'d_m',d_m,'cbeta',cbeta,...
                   'rprim',rprim,'mu_m',mu_m,'ilim_m',ilim_m,'Vlim_m',Vlim_m);
Deltat = struct('simcir',dt_simcir,'simvel',dt_simvel,'simpos',dt_simpos);
sys = ss(A,B,C,D);
dsys = c2d(sys,Deltat.simvel);
Matrix = struct('A',A,'B',B,'dA',dsys.a,'dB',dsys.b);
clear dsys sys sinalpha cosalpha alpha cbeta mB2 mD mR r C1 C2 ;
clear Dc dc etac I I0 ilimc Gc JAc Jc KEc Komega KTc Lc M m0 mc muc...
    numi rc Rc rgast thetagast ROT rprim Vlimc;
clear -regexp mK mI mM _m mA;



