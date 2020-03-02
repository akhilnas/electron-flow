classdef ClassSubstrate
    
   properties
      bias = {};
      Schottky = 0;
      Ohmic = 0;
      Zero_Slope = 0;
   end
   
   methods
       function r = CheckValues(Obj)
           r = 0;
           if ~Obj.Ohmic && ~Obj.Schottky && ~Obj.Zero_Slope
               msgStr = 'error: Could not find the barrier type / boundary condition for ''Substrate'' ';
               err    = MException('MATLAB:ValueCheckError',msgStr);
               throw(err);
           end
           r = 1;
       end
   end
   
end