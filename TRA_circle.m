classdef TRA_circle
   properties (SetAccess = public, GetAccess = public)
        R;
        alpha;
        phi0;
        x0;
        y0;
        theta0;
   end
   properties(SetAccess = private, GetAccess = public)
   end
   methods
       function obj =TRA_circle(R,x0,y0,theta0,phi,alpha)
           obj.R=R;
           obj.x0=x0;
           obj.y0=y0;
           obj.phi0 = phi;
           obj.theta0 = theta0;
           obj.alpha = alpha;
       end
       function obj = setorigin(obj,x,y,theta,phi)
           obj.x0 = x + obj.R*cos(phi);
           obj.y0 = y + obj.R*sin(phi);
           obj.phi0 = (pi+phi);
           obj.theta0 = theta;
       end
       function p = p(obj,s)
           x = obj.x0 + obj.R*cos(s+obj.phi0);
           y = obj.y0 + obj.R*sin(s+obj.phi0);
           theta = obj.theta0 + obj.alpha*s;
           p = [x;y;theta];
       end
       function v = v(obj,s,vs)
           vx = -obj.R*vs*sin(s+obj.phi0);
           vy = obj.R*vs*cos(s+obj.phi0);
           vtheta = obj.alpha*vs;
           v = [vx;vy;vtheta];
       end
   end
end