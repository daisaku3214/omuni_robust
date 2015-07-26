%%
%robot_define
%this file is definition of robot parameter
%%
%common definition
dt_simcir =1e-4;    %sampling time for circuit equation
dt_simvel =1e-3;    %sampling time for velocity simulation
dt_simpos =2.5e-2;    %sampling time for position simulation
g = 9.80665;        %[m/s^2] the Acceleration of gravity
ell = 4;            %the number of motor
I0 = 6;             %[kgm^2]    z-axis inertia of robot other than motor units
m0 = 20;             %[kg]       mass of robot other than motor units
rgast = 0.25;       %[m]        the radius between COM of robot and robot origin
thetagast = pi/4;   %[rad]      the radian from x ast axis
alpha = linspace(0,2*pi-2*pi/ell,ell);
rc = 0.5;           %[m]        the radius between omni wheel and robot origin
r = rc*ones(1,ell);
rprim = sqrt(rgast^2.+r.^2 -2*rgast.*r.*cos(thetagast-alpha));
cbeta = (rprim./r + r./rprim -rgast.^2./(r.*rprim))/2;
%sensor definition
num_i = 0;          %the number of current sensor
num_mring = 2;      %the number of measurement rings
num_jairo = 1;      %the number of jairo sensor
num_taco = ell;     %the number of taco meter by rotary encorder
num_LRF = 1;        %the number of LRF sensor
idynac = 4096;      %[bit] the bit number of current sensor
irangec = [-10 10]; %[A]   the dynamic range of current sensor
mringDc = 0.08;     %[m]   the diameter of measurement rings
mringc=2*pi/(4*100)...
 *ones(1,num_mring);%[rad/bit] the convert constants rad to bit
mring_alpha=[0 pi/2];%[rad] the place of measurment rings
mringrc = [0.1 0.1];%[m] the place of measurment rings
mring_dir = [0 pi/2];%[rad] the direction of mesurement rings
jairodynac=1024;    %[bit] the bit number of jairo sensor
jairorangec=[-2 2];%[rad/sec] the dynamic range of jairo sensor
jairoc = (jairorangec(2)-jairorangec(1))/jairodynac;
                    %[rad/sec.bit] the convert constans rad/s to bit
jairozero = jairodynac/2;%[bit]the zero point at jairo
tacoc=2*pi/(4*500)...
  *ones(1,num_taco);%[bit/rad] the convet constants rad to bit
LRF_dt=dt_simpos;   %[sec] the sampling time of LRF
LRFnoise = 0.0001*diag([1 1 pi]);%the Variance-covariance matrix of LRF noise
ENCnoise = 0.000001*diag([1 1 pi]);%the Variance-covariance matrix of encorder noise
TACnoise = 0.000001*diag([1 1 pi]); %the Variance-covariance matrix of tacometer noise
JROnoise = 0.00001; %the covariansce of jairo noise
JROoff = 1e-9;%[rad] the offset noise ratio of jairo
%%
%motor unit definition
%at this mode, all motor is same
%this patteren, I will use RE40 with 203113(1:3.5,eta = 90%)
Rc = 0.299;         %[ohm]      The internal resistance of the nominal motor
Lc = 0.0823/1000;   %[H]        The internal reactance of the nominal motor
KTc = 30.2/1000;    %[Nm/A]     The Torque constant of the nominal motor
ROT = 317;          %[rpm/V]    The Rotational constant of the nominal motor
ROT=2*pi*ROT/60;    %[rad/Vs]   The Rotational constant of the nominal motor
KEc = 1/(ROT);      %[Vs/rad]   The Back electromotive force constanto of the nominal motor
JAc = 142;          %[gcm^2]    The inertia of the rotor of the nominal motor
JAc =JAc/10^(7);    %[kgm^2]    The inertia of the rotor of the nominal motor
JGhc = 20+9.1;      %[gcm^2]    The inertia of the gear-head
JGhc = JGhc/10^(7); %[kgm^2]    The inertia of the gear-head
Gghc = 4.8;          %           The ratio of gear-head
Gcupc = 2.1;        %           The ratio of cuppring gear
Gc =  Gghc*Gcupc;    %           The ratio of gear
etac = 0.9*0.8;     %           The efficiency of gear
muc = 0.4;          %           The Coefficient of static friction
ilimc = 10;         %[A]        the limit of motor current
Vlimc = 22;         %[V]        the limit of input voltage
%this pattern, I use omni wheel of Diameter 185mm
Dc = 0.170;         %[m]        The dirmeter of omuni wheel
hc = 0.045;         %[m]        The width of omuni wheel
Jc = 0.213;         %[kgm^2]    The inertia of the omni whell
mc = 1.5;           %[kg]       The mass of a motor unit
dc = 0.01;          %[Nm/sec]    The declease ratio

R_m = Rc*ones(1,ell);
L_m = Lc*ones(1,ell);
KT_m = KTc*ones(1,ell);
KE_m = KEc*ones(1,ell);
JA_m = JAc*ones(1,ell);
G_m = Gc*ones(1,ell);
eta_m = etac*ones(1,ell);
D_m = Dc*ones(1,ell);
h_m = hc*ones(1,ell);
J_m = (Jc+Gcupc^2*JGhc)*ones(1,ell);
m_m = mc*ones(1,ell);
d_m = dc*ones(1,ell);
mu_m = muc*ones(1,ell);
ilim_m=ilimc*ones(1,ell);
Vlim_m=Vlimc*ones(1,ell);

M = m0 + sum(m_m);
I = I0 + sum(rprim.^2.*m_m);
%%
%define of uncertain
unnum_m = 10;       %the number of uncertain parameters by motor
unnum_o = 1;        %the number of unsertain parameters by other motor
unnum = ell*unnum_m...
      + unnum_o;    %the number of uncertain parameters
ungainnom = 1;
ungainper = [-20 20];
for i=1:ell

end
unLRF = LRFnoise;%the uncertain conversion matrix of LRF
unENC = ENCnoise;%the uncertain conversion matrix of encoder
unTAC = TACnoise;%the uncertain conversion matrix of tacometer
unJRO = JROnoise;%the uncertain of jairo
unJROoff=JROoff; %the uncertain of jairo offset
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
mAi = [-mMAs;mMAc;mIAr];
mAA = mM\mAi;
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
C2 = [zeros(num_i,3) eye(num_i,ell)];
C = [C1;C2];
D = zeros(size(C,1),ell);
v2omega = C1(2:ell+1,1:3);
omega2v = pinv(v2omega);
omega2Vol = diag(KE_m.*G_m);


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


%%
parameter = struct('g',g,'ell',ell,'I0',I0,'m0',m0,'rgast',rgast,...
                   'thetagast',thetagast,'alpha',alpha,'r',r,...
                   'num_i',num_i,'R_m',R_m,'L_m',L_m,'KT_m',KT_m,'KE_m',KE_m,...
                   'JA_m',JA_m,'G_m',G_m,'eta_m',eta_m,'D_m',D_m,'h_m',h_m,...
                   'J_m',J_m,'m_m',m_m,'d_m',d_m,'cbeta',cbeta,...
                   'rprim',rprim,'mu_m',mu_m,'ilim_m',ilim_m,'Vlim_m',Vlim_m);
Deltat = struct('simcir',dt_simcir,'simvel',dt_simvel,'simpos',dt_simpos);
sensparas = struct('num_i',num_i,'num_mring',num_mring,'num_jairo',num_jairo,...
          'num_taco',num_taco,'num_LRF',num_LRF,'idynac',idynac,'irangec',irangec,...
          'mringDc',mringDc,'mringc',mringc,'mring_alpha',mring_alpha,...
          'mringrc',mringrc,'mring_dir',mring_dir,'jairodynac',jairodynac,'jairozero',jairozero,...
          'jairorangec',jairorangec,'jairoc',jairoc,'tacoc',tacoc,'LRF_dt',LRF_dt);
Uncertain = struct('unnum',unnum,'unLRF',unLRF,'unENC',unENC,'unJRO',unJRO,...
           'unJROoff',unJROoff,'unTAC',unTAC);
sys = ss(A,B,C,D);
dsys = c2d(sys,Deltat.simvel);
Matrix = struct('A',A,'B',B,'dA',dsys.a,'dB',dsys.b);
clear dsys sys sinalpha cosalpha alpha cbeta mB2 mD mR r C1 C2 ;
clear Dc dc etac I I0 ilimc Gc JAc Jc KEc Komega KTc Lc M m0 mc muc...
    numi rc Rc rgast thetagast ROT rprim Vlimc;
clear -regexp mK mI mM _m mA;



