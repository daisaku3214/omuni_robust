%%
%robot_define
%this file is definition of robot parameter
%%
%common definition
dt_simcir =1e-5;%sampling time for circuit equation
dt_simvel =1e-3;%sampling time for velocity simulation
dt_simpos =1e-2;%sampling time for position simulation
ell = 4;        %the number of motor
I0 = 10;         %[kgm^2]    z-axis inertia of robot other than motor units
m0 = 10;        %[kg]       mass of robot other than motor units
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
Rc = 0.314;     %[ohm]      The internal resistance of the nominal motor
Lc = 0.085/1000;%[H]        The internal reactance of the nominal motor
KTc = 19.4/1000;%[Nm/A]     The Torque constant of the nominal motor
ROT = 491;      %[rpm/V]    The Rotational constant of the nominal motor
ROT=2*pi*ROT/60;%[rad/Vs]   The Rotational constant of the nominal motor
KEc = 1/(ROT);  %[Vs/rad]   The Back electromotive force constanto of the nominal motor
JAc = 68.1;     %[gcm^2]    The inertia of the rotor of the nominal motor
JAc =JAc/10^(7);%[kgm^2]    The inertia of the rotor of the nominal motor
Gc =  50;       %           The ratio of gear
etac = 0.9;     %           The efficiency of gear
%this pattern, I use omni wheel of Diameter 200mm
Dc = 0.3;       %[m]        The dirmeter of omni wheel
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

M = m0 + sum(m_m);
I = I0 + sum(rprim.^2.*m_m);
%%
sinalpha = sin(alpha);
cosalpha = cos(alpha);
a_m = 2*eta_m.*G_m.*KT_m./D_m;
b_m = 2*(J_m+eta_m.*G_m.^2.*JA_m)./D_m;
c_m = 2*b_m./D_m;
mMCss = sum(c_m.*sinalpha.^2)/M;
mMCcc = sum(c_m.*cosalpha.^2)/M;
mMCsc = sum(c_m.*sinalpha.*cosalpha)/M;
mMCrs = sum(r.*c_m.*sinalpha)/M;
mMCrc = sum(r.*c_m.*cosalpha)/M;
mICbs = sum(rprim.*c_m.*cbeta.*sinalpha)/I;
mICbc = sum(rprim.*c_m.*cbeta.*cosalpha)/I;
mICbr = sum(rprim.*r.*c_m.*cbeta)/I;
mMAs = a_m.*sinalpha/M;
mMAc = a_m.*cosalpha/M;
mIAr = rprim.*a_m.*cbeta/I;
mM = [1-mMCss mMCsc mMCrs;
      mMCsc 1-mMCcc -mMCrc;
      mICbs -mICbc 1-mICbr];
mK = [-mMAs;mMAc;mIAr];
mA1 = mM\mK;
mR  = diag(R_m./L_m);
wR  = mR;
mB2 = diag(ones(1,ell)./L_m);
wB2 = mB2;
Komega = 2*G_m.*KE_m./(L_m.*D_m);
mKx = Komega.*sin(alpha);
mKy = Komega.*cos(alpha);
mKq = Komega.*r;
mComega = [-mKx' mKy' mKq'];
mA = [zeros(3,3) mA1;
      mComega -mR];
mB = [zeros(3,ell);mB2];

mC1 = [0 0 1 zeros(1,ell);
       (2*sinalpha./D_m)' -(2*cosalpha./D_m)' -(2*rc./D_m)' zeros(ell,ell)];
mC2 = [zeros(numi,3) eye(numi,ell)];
mC = [mC1;mC2];
mD = zeros(size(mC,1),ell);
U2omega = mC1(2:ell+1,1:3);
omega2U = pinv(U2omega);

disp(size(mA));
[~,~,~,~,k] = ctrbf(mA,mB,mC);
disp(sum(k));
[~,~,~,~,k] = obsvf(mA,mB,mC);
disp(sum(k));
lambda = eig(mA);
[~,indextemp] = sort(real(lambda),'descend');
lambda = lambda(indextemp);
clear indextemp;
lambdatemp = 0;
for i = 1:size(mA,1);
if (lambda(i) > 0)
    lambdatemp = lambdatemp + lambda(i);
end
end
lambda = lambdatemp;
clear lambdatemp;
disp(rank([mC;mA-lambda]));


parameter = struct('ell',ell,'I0',I0,'m0',m0,'rgast',rgast,...
                   'thetagast',thetagast,'alpha',alpha,'r',r,...
                   'numi',numi,'R_m',R_m,'L_m',L_m,'KT_m',KT_m,...
                   'KE_m',KE_m,'JA_m',JA_m,'G_m',G_m,'eta_m',eta_m,...
                   'D_m',D_m,'J_m',J_m,'m_m',m_m,'d_m',d_m);
Deltat = struct('simcir',dt_simcir,'simvel',dt_simvel,'simpos',dt_simpos);
