classdef omuni_controller
   properties (SetAccess = public, GetAccess = public)
       const;   %constant of calcuated by parameter
       para;    %parameter of this omuni robot
       sensparas;%parameter os this omuni robot's sensor
       sysMat;  %system matrix
       contpara;
       hatx;
       x;
       u;
   end
   properties(SetAccess = private, GetAccess = public)
       Deltat;
       pLRF;
       dpsixy;
       Kvd;
       Kpd;
       Kd;
   end
   methods
       function obj = omuni_controller(parameter,Deltat,Matrix,sensparas,controller_parameter,hatx0,u0)
          obj.Deltat = Deltat;
          obj.sysMat = Matrix;
          obj.para = parameter;
          Cmring = 0.5*[(cos(sensparas.mring_dir)./sensparas.mringDc)'...
                         (sin(sensparas.mring_dir)./sensparas.mringDc)'...
                         (sensparas.mringrc.*...
                         sin(sensparas.mring_dir-sensparas.mring_alpha)...
                         ./sensparas.mringDc)'];
          Ctaco = 2*[(obj.para.G_m.*cos(parameter.alpha)./parameter.D_m)'...
                      -(obj.para.G_m.*sin(parameter.alpha)./parameter.D_m)'...
                      -(obj.para.G_m.*parameter.r./parameter.D_m)'];           
          obj.const = struct('Cmring',Cmring,'invCmring',pinv(Cmring),'Ctaco',Ctaco,'invCtaco',pinv(Ctaco));

          obj.hatx = hatx0;
          obj.pLRF = hatx0(1:3)*ones(1,sensparas.num_LRF);
          obj.dpsixy = zeros(2,1);
          obj.x  = hatx0;
          obj.u = u0;
          obj.sensparas =sensparas;
          obj.contpara = controller_parameter;
       end
       function obj = setKvdKpd(obj,Kvd,Kpd)
           obj.Kvd = Kvd;
           obj.Kpd = Kpd;
       end
       function obj = setKd(obj,Kd)
           obj.Kd = Kd;
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
       function T = Ttranshat(obj)
           sinq = sin(obj.hatx(3));
           cosq = cos(obj.hatx(3));
           T = [cosq sinq 0;
               -sinq cosq 0;
                  0    0  1];
       end
       function invT = invTtranshat(obj)
           sinq = sin(obj.hatx(3));
           cosq = cos(obj.hatx(3));
           invT = [cosq -sinq 0;
                   sinq  cosq 0;
                     0     0  1];
       end
       
       function hatvw = v_rev(obj,vw,dt)
           hatvw = zeros(size(vw));
           det = 1+(vw(3,1)*dt^2/2)^2;
           err = vw(3,1)*dt^2/2;

           hatvw(1,1) = vw(1,1)/det -err*vw(2,1)/det;
           hatvw(2,1) = err*vw(1,1)/det + vw(2,1)/det;
           hatvw(3,1) = vw(3,1);
       end
       function obj = calc_hatx(obj,sensdata)
           invT = obj.invTtranshat;
           newx = zeros(6,1);
           psi = invT*obj.const.invCmring*((obj.sensparas.mringc)'.*sensdata.dpsi);
           phi = invT*obj.const.invCtaco*((obj.sensparas.tacoc)'.*sensdata.dphi);
           omega = obj.sensparas.jairoc.*(sensdata.omega-obj.sensparas.jairozero);
           newx(6) = obj.contpara.fro*obj.Deltat.simvel*...
                   ((1-obj.contpara.alfENCvsJRO)*omega+obj.contpara.alfENCvsJRO*psi(3))...
                   +(1-obj.contpara.fro*obj.Deltat.simvel)*obj.hatx(6);

           newx(3) = obj.hatx(3) + newx(6)*obj.Deltat.simvel;
           obj.pLRF(3) = obj.pLRF(3)+newx(6)*obj.Deltat.simvel;
           
           obj.dpsixy = obj.contpara.fpxy*obj.Deltat.simvel* psi(1:2)+(1-obj.contpara.fpxy*obj.Deltat.simvel)*obj.dpsixy;

           newx(1:2) = obj.hatx(1:2)+ obj.dpsixy;
           obj.pLRF(1:2) = obj.pLRF(1:2) + obj.dpsixy;
           
%            newx(1:2) = obj.hatx(1:2)+ psi(1:2);
%            obj.pLRF(1:2) = obj.pLRF(1:2) + psi(1:2);
           
           if(sensdata.LRFflag==1)
             obj.pLRF = sensdata.pLRF;
%            newx(1:3) = (obj.contpara.alfp)*obj.pLRF+(1-(obj.contpara.alfp))*newx(1:3);
           end
           newx(1:3) = (obj.contpara.alfp/obj.contpara.LRFratio)*obj.pLRF ...
                     +(1-(obj.contpara.alfp/obj.contpara.LRFratio))*newx(1:3);
           newx(4:5) =obj.contpara.fvxy*phi(1:2)+(1-obj.Deltat.simvel*obj.contpara.fvxy)*obj.hatx(4:5);
           newi = obj.sysMat.dA(4:end,:)*[newx(4:6);obj.hatx(7:end)]...
                + obj.sysMat.dB(4:end,:)*obj.u;
           obj.hatx = [newx;newi];
       end
       function obj = Cont_LQR_cas(obj,xref,u_ff)
           Ttrans = obj.Ttranshat;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.hatx(1:3) - xref(1:3));
           ve =  Ttrans*obj.hatx(4:6) - vast;
           vref = [vast;zeros(obj.para.ell,1)] + obj.Kpd*([xe;ve;obj.hatx(7:end,1)]);
           vref = [obj.v_rev(vref(1:3,1),obj.Deltat.simvel);vref(4:end)];

           obj.u = -obj.Kvd*([Ttrans*obj.hatx(4:6);obj.hatx(7:end)]-vref)+u_ff;
       end
       function obj = Cont_LQR_cas_noest(obj,x,xref,u_ff)
           obj.x = x;
           Ttrans = obj.Ttrans;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.x(1:3) - xref(1:3));
           ve =  Ttrans*obj.x(4:6) - vast;
           vref = [vast;zeros(obj.para.ell,1)] + obj.Kpd*([xe;ve;obj.x(7:end,1)]);
           vref = [obj.v_rev(vref(1:3),obj.Deltat.simvel);vref(4:end)];

           obj.u = -obj.Kvd*([Ttrans*obj.x(4:6);obj.x(7:end)]-vref)+u_ff;
       end
       function obj = Cont_LQR_cas_norev(obj,xref,u_ff)
           Ttrans = obj.Ttranshat;
           vast = Ttrans*xref(4:6);
           xe =  Ttrans*(obj.hatx(1:3) - xref(1:3));
           ve =  Ttrans*obj.hatx(4:6) - vast;
           vref = [vast;zeros(obj.para.ell,1)] + obj.Kpd*([xe;ve;obj.hatx(7:end,1)]);
           
           obj.u = -obj.Kvd*([Ttrans*obj.hatx(4:6);obj.hatx(7:end)]-vref)+u_ff;
       end
       function obj = Cont_LQR_cas_noest_norev(obj,x,xref,u_ff)
           obj.x = x;
           Ttrans = obj.Ttrans;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.x(1:3) - xref(1:3));
           ve =  Ttrans*obj.x(4:6) - vast;
           vref = [vast;zeros(obj.para.ell,1)] + obj.Kpd*([xe;ve;obj.x(7:end,1)]);
           
           obj.u = -obj.Kvd*([Ttrans*obj.x(4:6);obj.x(7:end)]-vref)+u_ff;
       end
       function obj = Cont_LQR_cas_filt(obj,xref,u_ff)
           Ttrans = obj.Ttranshat;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.hatx(1:3) - xref(1:3));
           ve =  Ttrans*obj.hatx(4:6) - vast;
           vref = [vast;zeros(obj.para.ell,1)] + obj.Kpd*([xe;ve;obj.hatx(7:end,1)]);
           vref = [obj.v_rev(vref(1:3),obj.Deltat.simvel);vref(4:end)];

           utemp = -obj.Kvd*([Ttrans*obj.hatx(4:6);obj.hatx(7:end)]-vref);
           obj.u = obj.Deltat.simvel*obj.contpara.fu*utemp...
                   +(1-obj.Deltat.simvel*obj.contpara.fu)*(obj.u-u_ff)+u_ff;
       end
       function obj = Cont_LQR_cas_norev_filt(obj,xref,u_ff)
           Ttrans = obj.Ttranshat;
           vast = Ttrans*xref(4:6);
           xe =  Ttrans*(obj.hatx(1:3) - xref(1:3));
           ve =  Ttrans*obj.hatx(4:6) - vast;
           vref = [vast;zeros(obj.para.ell,1)] + obj.Kpd*([xe;ve;obj.hatx(7:end,1)]);
           
           utemp = -obj.Kvd*([Ttrans*obj.hatx(4:6);obj.hatx(7:end)]-vref);
           obj.u = obj.Deltat.simvel*obj.contpara.fu*utemp...
                   +(1-obj.Deltat.simvel*obj.contpara.fu)*(obj.u-u_ff)+u_ff;
       end

       
       function obj = Cont_LQR_single(obj,xref,u_ff)
           Ttrans = obj.Ttranshat;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.hatx(1:3) - xref(1:3));
           ve =  Ttrans*obj.hatx(4:6) - vast;
           obj.u = -obj.Kd*([xe;ve;obj.hatx(7:end)])+u_ff;
       end
       function obj = Cont_LQR_single_filt(obj,xref,u_ff)
           Ttrans = obj.Ttranshat;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.hatx(1:3) - xref(1:3));
           ve =  Ttrans*obj.hatx(4:6) - vast;
           utemp = -obj.Kd*([xe;ve;obj.hatx(7:end)]);
           obj.u = obj.Deltat.simvel*obj.contpara.fu*utemp...
                   +(1-obj.Deltat.simvel*obj.contpara.fu)*(obj.u-u_ff)+u_ff;
       end
       function obj = Cont_LQR_single_noest(obj,x,xref,u_ff)
           obj.x = x;
           Ttrans = obj.Ttrans;
           vast = Ttrans*xref(4:6);
           xe =   Ttrans*(obj.x(1:3) - xref(1:3));
           ve =  Ttrans*obj.x(4:6) - vast;
           obj.u = -obj.Kd*([xe;ve;obj.x(7:end)])+u_ff;
       end       
   end
end