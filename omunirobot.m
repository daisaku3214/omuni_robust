classdef omunirobot
   properties (SetAccess = public, GetAccess = public)
       para;    %parameter of this omuni robot
       const;   %constant of calcuated by parameter
       sensparas;%parameter os this omuni robot's sensor
       sysMat;  %system matrix
       x;       %state variable of this omuni robot
                %x = [xy\theta\dot{x}\dot{y}\dot{\theta}i_1\cdots i_\ell]^T
                %size(x) = [6+para.ell 1];
       u;       %input variable of this omuni robot
                %u = [V_1\cdots V_\ell]^T
                %size(u) = [para.ell 1];
       F;       %grip force
       t;       %the time
       psi;     %[rad] the angles of mesure rings
       phi;     %[rad] the angles of motor meter
       dpsi;    %[bit] the angles of mesure rings by enc at each sampling time
       dphi;    %[bit] the angles of motor meter by enc at each sampling time
       pLRF;    %[m/rad] the posture parameter by LRF
       LRFflag; %[bool] flag of wheter or not LRF data comming
       omega;   %[rad] the anguler of rotation velosity at each sampling time
   end
   properties(SetAccess = private, GetAccess = public)
       LRFcounter;
       Deltat;
       Uncertain;
       sensdata;
   end
   methods
       function obj = omunirobot(parameter,Deltat,Matrix,Uncertain,sensparas,x0,u0)
           obj.para = parameter;
           obj.Deltat = Deltat;
           obj.sysMat = Matrix;
           obj.Uncertain =Uncertain;
           obj.sensparas =sensparas;
           a_m = 2*parameter.eta_m.*parameter.G_m.*parameter.KT_m./parameter.D_m;
           b_m = 2*(parameter.J_m+parameter.eta_m.*parameter.G_m.^2.*parameter.JA_m)./parameter.D_m;
           c_m = 2*b_m./parameter.D_m;
           dD_m =2*parameter.d_m./parameter.D_m;
           Komega_m = (2*parameter.G_m.*parameter.KE_m)./(parameter.L_m.*parameter.D_m);
           Ktheta_m = Komega_m.*parameter.r;
           wR = diag(parameter.R_m./parameter.L_m);
           wB2 = diag(1./parameter.L_m);
           wB = [zeros(6,parameter.ell);wB2];
           ratiocir = Deltat.simvel/Deltat.simcir;
           ratiovel = Deltat.simpos/Deltat.simvel;
           Hinv = pinv([ones(1,parameter.ell);
                 parameter.r.*sin(parameter.alpha);
                 parameter.r.*cos(parameter.alpha)]);
           N = Hinv(:,1)*(parameter.m0+sum(parameter.m_m))*parameter.g;
           ilim = min(abs(parameter.ilim_m));
           Vlim = min(abs(parameter.Vlim_m));
           ratioLRF = sensparas.LRF_dt/Deltat.simvel;

           Cmring = 0.5*[(cos(sensparas.mring_dir)./sensparas.mringDc)'...
                         (sin(sensparas.mring_dir)./sensparas.mringDc)'...
                         (sensparas.mringrc.*...
                         sin(sensparas.mring_dir-sensparas.mring_alpha)...
                         ./sensparas.mringDc)'];
           Ctaco = 2*[(obj.para.G_m.*cos(parameter.alpha)./parameter.D_m)'...
                      -(obj.para.G_m.*sin(parameter.alpha)./parameter.D_m)'...
                      -(obj.para.G_m.*parameter.r./parameter.D_m)'];
           obj.const = struct('a_m',a_m,'b_m',b_m,'c_m',c_m,'dD_m',dD_m,...
                              'Komega_m',Komega_m,'Ktheta_m',Ktheta_m,...
                              'wR',wR,'wB2',wB2,'wB',wB,'ratiocir',ratiocir,...
                              'ratiovel',ratiovel,'N',N,'ilim',ilim,...
                              'Vlim',Vlim,'ratioLRF',ratioLRF,'Cmring',Cmring,...
                              'invCmring',pinv(Cmring),'Ctaco',Ctaco,'invCtaco',pinv(Ctaco));
           obj.x = x0;
           obj.u = u0;
           obj.t = 0;
           obj.F = obj.calcF(zeros(obj.para.ell));

           obj.LRFcounter = 1;
           
           obj.psi = zeros(sensparas.num_mring,1);
           obj.phi = zeros(sensparas.num_taco,1);
           
           obj.dpsi = zeros(sensparas.num_mring,1);
           obj.dphi = zeros(sensparas.num_taco,1);
           obj.pLRF = x0;
           obj.LRFflag = 1;
           obj.omega = zeros(sensparas.num_jairo,1);
           obj.sensdata = struct('dpsi',obj.dpsi,'dphi',obj.dphi,...
                                 'pLRF',obj.pLRF,'LRFflag',obj.LRFflag,...
                                 'omega',obj.omega);
       end

       function T = Ttrans(obj)
           sinq = sin(obj.x(3));
           cosq = cos(obj.x(3));
           T = [cosq sinq 0;
               -sinq cosq 0;
                  0    0  1];
       end
       function invT = invTtrans(obj)
           sinq = sin(obj.x(3));
           cosq = cos(obj.x(3));
           invT = [cosq -sinq 0;
                   sinq  cosq 0;
                     0     0  1];
       end
       function Tarray = Tarray(obj,theta,inputs)
           sinq = sin(theta(1,:));
           cosq = cos(theta(1,:));
           length = size(theta,2);
           Tarray = zeros(3,3);
           for i = 1:length
              Trans =  [cosq(i) sinq(i) 0;
                        -sinq(i)  cosq(i) 0;
                          0        0     1]; 
              Tarray(:,i) = Trans*inputs(1:3,i);
           end
           
       end
       function limitu = Vlimitation(obj,u)
           if (max(abs(u)) > obj.const.Vlim)
               limitu = obj.const.Vlim*u/max(abs(u));
           else
               limitu = u;
           end
       end
       function omega = calcomega(obj)
           theta_w = obj.x(3) + obj.para.alpha;
           sinqw = sin(theta_w);
           cosqw = cos(theta_w);
           o_x = 2*sinqw./obj.para.D_m;
           o_y = 2*cosqw./obj.para.D_m;
           
           omega = [o_x' -o_y' -obj.para.r']*obj.x(4:6);
       end
       function [dx,F] = calcdx(obj)
           theta_w = obj.x(3) + obj.para.alpha;
           sinqw = sin(theta_w);
           cosqw = cos(theta_w);
           sinqg = sin(obj.x(3)+obj.para.thetagast);
           cosqg = cos(obj.x(3)+obj.para.thetagast);
           rprim = sqrt(obj.para.rgast^2 + obj.para.r.^2 ...
           -2*obj.para.rgast.*obj.para.r.^2.*cos(obj.para.thetagast - obj.para.alpha));
           cbeta = rprim./(2*obj.para.r) +obj.para.r./(2*rprim) ...
               -obj.para.rgast./(2*rprim.*obj.para.r);
           I = obj.para.I0 + sum(obj.para.m_m.*rprim.^2);
           M = obj.para.m0 + sum(obj.para.m_m);

           ssqw = sinqw.^2;
           scqw = sinqw.*cosqw;
           ccqw = cosqw.^2;
           cbsqw = cbeta.*sinqw;
           cbcqw = cbeta.*cosqw;
           
           Kx_m = obj.const.Komega_m.*sinqw;
           Ky_m = obj.const.Komega_m.*cosqw;
           
           wMCss = sum(obj.const.c_m.*ssqw)/M;
           wMCcc = sum(obj.const.c_m.*ccqw)/M;
           wMCsc = sum(obj.const.c_m.*scqw)/M;
           wMCrs = sum(obj.para.r.*obj.const.c_m.*sinqw)/M;
           wMCrc = sum(obj.para.r.*obj.const.c_m.*cosqw)/M;
           wICbs = sum(rprim.*obj.const.c_m.*cbsqw)/I;
           wICbc = sum(rprim.*obj.const.c_m.*cbcqw)/I;
           wICbr = sum(rprim.*obj.para.r.*obj.const.c_m.*cbeta)/I;
           wMDss = sum(obj.const.dD_m.*ssqw)/M;
           wMDcc = sum(obj.const.dD_m.*ccqw)/M;
           wMDsc = sum(obj.const.dD_m.*scqw)/M;
           wMDrs = sum(obj.para.r.*obj.const.dD_m.*sinqw)/M;
           wMDrc = sum(obj.para.r.*obj.const.dD_m.*cosqw)/M;
           wIDbs = sum(rprim.*obj.const.dD_m.*cbsqw)/I;
           wIDbc = sum(rprim.*obj.const.dD_m.*cbcqw)/I;
           wIDbr = sum(rprim.*obj.para.r.*obj.const.dD_m.*cbeta)/I;
           wMAs = (obj.const.a_m.*sinqw)/M;
           wMAc = (obj.const.a_m.*cosqw)/M;
           wIAr = (rprim.*obj.const.a_m.*cbeta)/I;
           wM = [1-wMCss wMCsc wMCrs-obj.para.rgast*sinqg;
                 wMCsc 1-wMCcc -wMCrc+obj.para.rgast*cosqg;
                 wICbs -wICbc 1-wICbr];
           wC = [ wMCsc  wMCss obj.para.rgast*cosqg;
                 -wMCcc  -wMCsc obj.para.rgast*sinqg;
                 -wICbc  -wICbs 0];
           wD = [ wMDss -wMDsc -wMDrs;
                 -wMDsc  wMDcc  wMDrc;
                 -wIDbs  wIDbc  wIDbr];
           wAi = [-wMAs;wMAc;wIAr];
           wF1 = wM\((wC*obj.x(6)+wD)*obj.x(4:6)+wAi*obj.x(7:end));
           wKomega = [-Kx_m' Ky_m' obj.const.Ktheta_m'];
           wF2 = wKomega*obj.x(4:6) - obj.const.wR*obj.x(7:end)+...
               obj.const.wB2*obj.u;
           dx = [obj.x(4:6);wF1;wF2];
           dxq = -dx(5)+obj.x(4)*obj.x(6);
           dyq = dx(4)+obj.x(5)*obj.x(6);
           Fai = (obj.const.a_m.*obj.x(7:end)')';
           Fcx = -(obj.const.c_m.*(sinqw*(dyq)+cosqw*(dxq)-obj.para.r*dx(6)))';
           Fdx = -((obj.para.d_m./obj.para.D_m).*(sinqw*obj.x(4)...
                -cosqw*obj.x(5)-obj.para.r*obj.x(6)))';
            F = Fai+Fcx+Fdx;
       end
       function F = calcF(obj,dx)
           sinqw = sin(obj.x(3)+obj.para.alpha);
           cosqw = cos(obj.x(3)+obj.para.alpha);
           dxq = -dx(5)+obj.x(4)*obj.x(6);
           dyq = dx(4)+obj.x(5)*obj.x(6);
           Fai = (obj.const.a_m.*obj.x(7:end)')';
           Fcx = -(obj.const.c_m.*(sinqw*(dyq)+cosqw*(dxq)-obj.para.r*dx(6)))';
           Fdx = -((obj.para.d_m./obj.para.D_m).*(sinqw*obj.x(4)...
                -cosqw*obj.x(5)-obj.para.r*obj.x(6)))';
            F = Fai+Fcx+Fdx;
       end
       function obj = shiftx(obj)
           obj.u = obj.Vlimitation(obj.u);
           for i = 1:obj.const.ratiocir
           tempx = obj.x;
               [k1,Fk1] = obj.calcdx;
               obj.x = tempx + k1*obj.Deltat.simcir/2;
               [k2,Fk2] = obj.calcdx;
               obj.x = tempx + k2*obj.Deltat.simcir/2;
               [k3,Fk3] = obj.calcdx;
               obj.x = tempx + k3*obj.Deltat.simcir;
               [k4,Fk4] = obj.calcdx;
               dx = (k1+2*k2+2*k3+k4)*obj.Deltat.simcir/6;
               obj.x = tempx + dx;
%                clc;disp([Fk1 Fk2 Fk3 Fk4]);
               Ttrans = obj.Ttrans;
               obj.F = (Fk1+2*Fk2+2*Fk3+Fk4)/6;
               obj.psi =obj.psi+obj.const.Cmring*Ttrans*dx(1:3);
               obj.phi =obj.phi+obj.const.Ctaco*Ttrans*dx(1:3);
           end
           obj.t = obj.t + obj.Deltat.simvel;
           obj.LRFcounter = obj.LRFcounter+1;
           
           if(obj.LRFcounter==obj.const.ratioLRF)
               obj = obj.getLRF;
               obj.LRFflag = 1;
               obj.LRFcounter = 1;
           else
               obj.LRFflag = 0;
           end
           obj = obj.getTAC;
           obj = obj.getENC;
           obj = obj.getJRO;
           obj.sensdata.dpsi = obj.dpsi;
           obj.sensdata.dphi = obj.dphi;
           obj.sensdata.pLRF = obj.pLRF;
           obj.sensdata.LRFflag = obj.LRFflag;
           obj.sensdata.omega = obj.omega;
       end
       
       function obj= getLRF(obj)
           obj.pLRF = obj.x(1:3)+obj.Uncertain.unLRF*randn(3,1);
       end
       function obj = getTAC(obj)
           obj.dphi = fix(obj.phi./(obj.sensparas.tacoc)');
           obj.phi = obj.phi - obj.dphi.*(obj.sensparas.tacoc)';
           obj.dphi = obj.dphi + obj.const.Ctaco*obj.Uncertain.unTAC*randn(3,1);
       end
       function obj = getENC(obj)
           obj.dpsi = fix(obj.psi./(obj.sensparas.mringc)');
           obj.psi = obj.psi - obj.dpsi.*(obj.sensparas.mringc)';%restcounter
           obj.dpsi = obj.dpsi+ obj.const.Cmring*obj.Uncertain.unENC*randn(3,1);
       end
       function obj = getJRO(obj)
           tempomega = obj.x(6) + obj.Uncertain.unJRO*randn(1)...
                 + obj.Uncertain.unJROoff*obj.t;
           if(tempomega < obj.sensparas.jairorangec(1))
              tempomega = obj.sensparas.jairorangec(1);
           end
           if(tempomega > obj.sensparas.jairorangec(2))
              tempomega = obj.sensparas.jairorangec(1);
           end
           obj.omega = fix(tempomega./obj.sensparas.jairoc)+obj.sensparas.jairozero;
       end

   end
end