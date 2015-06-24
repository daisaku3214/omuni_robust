classdef omuniplants
   properties (SetAccess = public, GetAccess = public)
       numP;
       nominal;
       real;
   end
      methods
       function obj = omuniplants(parameter,Deltat,x0,u0,numP)
           obj.numP = numP;
           obj.nominal = omunirobot(parameter,Deltat,x0,u0);
           obj.real = obj.nominal;
           for i = 1:numP
             newreal = obj.nominal;
             newreal = newreal.restlog(x0+randn(size(newreal.x)),u0);
             obj.real = [obj.real;newreal];
           end
       end
       function obj = setnewuplants(obj,newu)
           obj.nominal.u = newu;
           for i = 1:obj.numP
               obj.real(i).u = newu;
           end
       end
       function obj = shiftplants(obj)
           obj.nominal = obj.nominal.shiftx;
           for i = 1:obj.numP
               obj.real(i) = obj.real(i).shiftx;
           end
       end

    end
end