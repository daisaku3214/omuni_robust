classdef TRA_line
   properties (SetAccess = public, GetAccess = public)
        x0;
        y0;
        theta0;
        phi0;
        alpha;
   end
   properties(SetAccess = private, GetAccess = public)
   end
   methods
       function obj =TRA_line(x0,y0,theta0,phi0,alpha)
           obj.x0=x0;
           obj.y0=y0;
           obj.theta0 = theta0;
           obj.phi0 = phi0;
           obj.alpha = alpha;
       end
       function p = p(obj,s)
           x = obj.x0 + s*cos(obj.phi0);
           y = obj.y0 + s*sin(obj.phi0);
           theta = obj.theta0 + obj.alpha*s;
           p = [x;y;theta];
       end
       function v = v(obj,~,vs)
           vx = vs*cos(obj.phi0);
           vy = vs*sin(obj.phi0);
           vtheta = obj.alpha*vs;
           v = [vx;vy;vtheta];
       end
   end
end