classdef omunimachine
    properties (SetAccess = public, GetAccess = public)
        robot;
        controller;
    end
    properties(SetAccess = private, GetAccess = public)
       Xlog;
       hatXlog;
       Xreflog;
       Ulog;
       nominalUlog;
       Flog;
       Tlog;
    end
    methods
       function obj = omunimachine(parameter,Deltat,Matrix,Uncertain,...
                      sensparas,controller_parameter,x0,hatx0,u0)
           obj.robot = omunirobot(parameter,Deltat,Matrix,Uncertain,sensparas,x0,u0);
           obj.controller = omuni_controller(parameter,Deltat,Matrix,sensparas,controller_parameter,hatx0,u0);

           obj.Xlog = x0;
           obj.hatXlog = hatx0;
           obj.Xreflog = hatx0;
           obj.Ulog = u0;
           obj.nominalUlog = u0;
           obj.Flog = obj.robot.F;
           obj.Tlog = obj.robot.t;
       end
       function obj = restlog(obj,x0,u0)
           obj.robot.x = x0;
           obj.robot.u = u0;
           obj.robot.nominalu = u0;
           obj.robot.t = 0;
           obj.robot.F = obj.calcF(zeros(obj.para.ell));
           obj.controller.hatx = x0;
           obj.controller.u=u0;

           obj.Xlog = x0;
           obj.hatXlog = hatx0;
           obj.Xreflog = hatx0;
           obj.Ulog = u0;
           obj.nominalUlog = u0;
           obj.Flog = obj.robot.F;
           obj.Tlog = obj.robot.t;
       end
       function obj = setKvdKpd(obj,Kvd,Kpd)
           obj.controller = obj.controller.setKvdKpd(Kvd,Kpd);
       end
       function obj = setKd(obj,Kd)
           obj.controller = obj.controller.setKd(Kd);
       end
       %%
       function obj = control_shift(obj,xref,u_ff,mode)
           
           obj.robot = obj.robot.shiftx;
           
           obj.controller = obj.controller.calc_hatx(obj.robot.sensdata);
                      
           obj.Xlog = [obj.Xlog obj.robot.x];
           obj.hatXlog = [obj.hatXlog obj.controller.hatx];
           obj.Xreflog = [obj.Xreflog xref];
           obj.Ulog = [obj.Ulog obj.robot.u];
           obj.nominalUlog = [obj.nominalUlog obj.controller.u];
           obj.Flog = [obj.Flog obj.robot.F];
           obj.Tlog = [obj.Tlog obj.robot.t];
           
           switch mode
               case 'LQR_single'
                   obj.controller = obj.controller.Cont_LQR_single(xref,u_ff);
               case 'LQR_single_filt'
                   obj.controller = obj.controller.Cont_LQR_single_filt(xref,u_ff);
               case 'LQR_single_noest'
                   obj.controller = obj.controller.Cont_LQR_single_noest(obj.robot.x,xref,u_ff);
               case 'LQR_cas'
                   obj.controller = obj.controller.Cont_LQR_cas(xref,u_ff);
               case 'LQR_cas_noest'
                   obj.controller = obj.controller.Cont_LQR_cas_noest(obj.robot.x,xref,u_ff);
               case 'LQR_cas_norev'
                   obj.controller = obj.controller.Cont_LQR_cas_norev(xref,u_ff);
               case 'LQR_cas_noest_norev'
                   obj.controller = obj.controller.Cont_LQR_cas_noest_norev(obj.robot.x,xref,u_ff);
               case 'LQR_cas_filt'
                   obj.controller = obj.controller.Cont_LQR_cas_filt(xref,u_ff);
               case 'LQR_cas_norev_filt'
                   obj.controller = obj.controller.Cont_LQR_cas_norev_filt(xref,u_ff);
           end
           
           obj.robot.u = obj.controller.u;

       end
       %%
       %plot log function
       function plotGrip(obj,index)
           Fmax = obj.robot.para.mu_m'.*obj.robot.const.N;
           plot(obj.Tlog,obj.Flog(index,:),'-m','LineWidth',2);
           hold on; grid on;
           xlim([0 obj.Tlog(1,end)]);
           plot([obj.Tlog(1,1) obj.Tlog(1,end)],[ Fmax(index)  Fmax(index)],':g','LineWidth',2);
           plot([obj.Tlog(1,1) obj.Tlog(1,end)],[-Fmax(index) -Fmax(index)],':g','LineWidth',2);
           ylabel(strcat(strcat('$$F_',num2str(index)),'[N]$$'), 'interpreter', 'latex');
           hold off;
       end
       %this function plot Current log(Xlog(7:end,;))
       %index is the motor index which you want to check current log
       %estimate is boolvalue of wheter show or not estimate current value
       function plotCurrent(obj,index,estimate)
           plot(obj.Tlog,obj.Xlog(6+index,:),'-m','LineWidth',2);
           hold on; grid on;
           if(estimate ==1)
               plot(obj.Tlog,obj.hatXlog(6+index,:),'-r','LineWidth',2);
           end
           xlim([0 obj.Tlog(1,end)]);
           plot([obj.Tlog(1,1) obj.Tlog(1,end)],[ obj.robot.para.ilim_m(index)  obj.robot.para.ilim_m(index)],':g','LineWidth',2);
           plot([obj.Tlog(1,1) obj.Tlog(1,end)],[-obj.robot.para.ilim_m(index) -obj.robot.para.ilim_m(index)],':g','LineWidth',2);
           ylabel(strcat(strcat('$$i_',num2str(index)),'[A]$$'), 'interpreter', 'latex');
           if(estimate ==1)
               legend(strcat('i_',num2str(index)),strcat('hat i_',num2str(index)));
           else
               legend(strcat('i_',num2str(index)));
           end
           hold off;
       end
       function plotVoltage(obj,index,nominal)
           plot(obj.Tlog,obj.Ulog(index,:),'-m','LineWidth',2);
           hold on; grid on;
           if(nominal ==1)
               plot(obj.Tlog,obj.nominalUlog(index,:),'-y','LineWidth',2);
               legend(strcat('V_',num2str(index)),strcat('nominalV_',num2str(index)));
           else
               legend(strcat('V_',num2str(index)));
           end
           xlim([0 obj.Tlog(1,end)]);
           ylabel(strcat(strcat('$$V_',num2str(index)),'[V]$$'), 'interpreter', 'latex');
           hold off;
       end
       function plotPosture(obj,index,estimate)
           plot(obj.Tlog,obj.Xlog(index,:),'-m','LineWidth',2);%posture
           hold on; grid on;
           if(estimate ==1)
               plot(obj.Tlog,obj.hatXlog(index,:),'-r','LineWidth',2);%estimate posutre
           end
           plot(obj.Tlog,obj.Xreflog(index,:),'--c','LineWidth',2)%refference
           if(estimate ==1)
               legend('x','hatx','xref');
           else
               legend('x','xref');
           end
           xlim([0 obj.Tlog(1,end)]);
           if index==1
               ylabel('$$x {[m]}$$', 'interpreter', 'latex');
           end
           if index==2
               ylabel('$$y {[m]}$$', 'interpreter', 'latex');
           end
           if index==3
               ylabel('$$\theta {[rad]}$$', 'interpreter', 'latex');
           end
           hold off;
       end
       function plotVelocity(obj,index,estimate)
           plot(obj.Tlog,obj.Xlog(3+index,:),'-m','LineWidth',2);%Velocity
           hold on; grid on;
           if(estimate ==1)
               plot(obj.Tlog,obj.hatXlog(3+index,:),'-r','LineWidth',2);%estimate Velocity
           end
           plot(obj.Tlog,obj.Xreflog(3+index,:),'--c','LineWidth',2);%refference
           if(estimate ==1)
               legend('x','hatx','xref');
           else
               legend('x','xref');
           end
           xlim([0 obj.Tlog(1,end)]);
           if index==1
               ylabel('$$\dot x {[m/s]}$$', 'interpreter', 'latex');
           end
           if index==2
               ylabel('$$\dot y {[m/s]}$$', 'interpreter', 'latex');
           end
           if index==3
               ylabel('$$\dot \theta {[rad/s]}$$', 'interpreter', 'latex');
           end
           hold off;
       end
       function plotPostureerror(obj,index,estimate)
           plot(obj.Tlog,obj.Xlog(index,:)-obj.Xreflog(index,:),'-m','LineWidth',2);
           hold on; grid on;
           if(estimate ==1)
               plot(obj.Tlog,obj.hatXlog(index,:)-obj.Xreflog(index,:),'-r','LineWidth',2);
           end
           plot([obj.Tlog(1,1) obj.Tlog(1,end)],zeros(1,2),'--c','LineWidth',2)
           if(estimate ==1)
               legend('x','hatx','xref');
           else
               legend('x','xref');
           end
           xlim([0 obj.Tlog(1,end)]);
           if index==1
               ylabel('$$e_x {[m]}$$', 'interpreter', 'latex');
           end
           if index==2
               ylabel('$$e_y {[m]}$$', 'interpreter', 'latex');
           end
           if index==3
               ylabel('$$e_\theta {[rad]}$$', 'interpreter', 'latex');
           end
           hold off;
       end
       function plotVelocityerror(obj,index,estimate)
           plot(obj.Tlog,obj.Xlog(index+3,:)-obj.Xreflog(index+3,:),'-m','LineWidth',2);
           hold on; grid on;
           if(estimate ==1)
               plot(obj.Tlog,obj.hatXlog(index+3,:)-obj.Xreflog(index+3,:),'-r','LineWidth',2);
           end
           plot([obj.Tlog(1,1) obj.Tlog(1,end)],zeros(1,2),'--c','LineWidth',2)
           if(estimate ==1)
               legend('x','hatx','xref');
           else
               legend('x','xref');
           end
           xlim([0 obj.Tlog(1,end)]);
           if index==1
               ylabel('$$e_{\dot{x}} {[m]}$$', 'interpreter', 'latex');
           end
           if index==2
               ylabel('$$e_{\dot{y}} {[m]}$$', 'interpreter', 'latex');
           end
           if index==3
               ylabel('$$e_{\dot{\theta}} {[rad]}$$', 'interpreter', 'latex');
           end
           hold off;
       end
       function plotPostureestimateerror(obj,index)
           plot(obj.Tlog,obj.Xlog(index,:)-obj.hatXlog(index,:),'-m','LineWidth',2);
           grid on;
           legend('x-hatx');
           xlim([0 obj.Tlog(1,end)]);
           if index==1
               ylabel('$$e_x {[m]}$$', 'interpreter', 'latex');
           end
           if index==2
               ylabel('$$e_y {[m]}$$', 'interpreter', 'latex');
           end
           if index==3
               ylabel('$$e_\theta {[rad]}$$', 'interpreter', 'latex');
           end
           hold off;
       end
       function plotVelocityestimateerror(obj,index)
           plot(obj.Tlog,obj.Xlog(index+3,:)-obj.hatXlog(index+3,:),'-m','LineWidth',2);
           grid on;
           legend('x-hatx');
           xlim([0 obj.Tlog(1,end)]);
           if index==1
               ylabel('$$e_x {[m]}$$', 'interpreter', 'latex');
           end
           if index==2
               ylabel('$$e_y {[m]}$$', 'interpreter', 'latex');
           end
           if index==3
               ylabel('$$e_\theta {[rad]}$$', 'interpreter', 'latex');
           end
           hold off;
       end
       
       function plotinputlogger(obj,estimate,nominal)
           %open figure window
           hfin = figure;
           set(hfin,'position',[150 150 1400 obj.robot.para.ell*150]);
           for i = 1:obj.robot.para.ell
               %plot grip force    
               subplot(obj.robot.para.ell,3,1+3*(i-1));
               obj.plotGrip(i);
               %plot motor current
               subplot(obj.robot.para.ell,3,2+3*(i-1));
               obj.plotCurrent(i,estimate);
               %plot input voltage
               subplot(obj.robot.para.ell,3,3+3*(i-1));
               obj.plotVoltage(i,nominal);
           end
           %set xlabel
           for j = 1:3
               subplot(obj.robot.para.ell,3,2*obj.robot.para.ell+1+j);
               xlabel('time [sec]', 'interpreter', 'latex');
           end
           %set graff title
           subplot(obj.robot.para.ell,3,2);
           title('Motor current', 'interpreter', 'latex');
           subplot(obj.robot.para.ell,3,2);
           title('Motor current', 'interpreter', 'latex');
           subplot(obj.robot.para.ell,3,3);
           title('Input voltage', 'interpreter', 'latex');           
       end
       function plotoutputlogger(obj,estimate)
           %open figure window
           hfout = figure;
           set(hfout,'position',[150 150 1400 600]);
           for i = 1:3
               %plot robot posture
               subplot(3,3,1+3*(i-1));
               obj.plotPosture(i,estimate);
               %plot robot velosity at world
               subplot(3,3,2+3*(i-1));
               obj.plotVelocity(i,estimate);
               subplot(3,3,3*i);
               %plot robot posture error
               obj.plotPostureerror(i,estimate);
           end
           %set xlabel
           for j = 1:3
               subplot(3,3,6+j);
               xlabel('time [sec]', 'interpreter', 'latex');
           end
           %set graff title
           subplot(3,3,1);
           title('Robot posture', 'interpreter', 'latex');
           subplot(3,3,2);
           title('Robot velosity at world', 'interpreter', 'latex');
           subplot(3,3,3);
           title('Robot posture error', 'interpreter', 'latex');

       end
       function plotestimateerrorlogger(obj)
           %open figure window
           hfestimate = figure;
           set(hfestimate,'position',[150 150 800 600]);
           for i = 1:3
               %plot robot posture
               subplot(3,2,1+2*(i-1));
               obj.plotPostureestimateerror(i);
               %plot robot velosity at world
               subplot(3,2,2+2*(i-1));
               obj.plotVelocityestimateerror(i);
           end
           %set xlabel
           for j = 1:2
               subplot(3,2,4+j);
               xlabel('time [sec]', 'interpreter', 'latex');
           end
           %set graff title
           subplot(3,2,1);
           title('Robot posture estimate error', 'interpreter', 'latex');
           subplot(3,2,2);
           title('Robot velosity estimate error at world', 'interpreter', 'latex');
       end
    end
end