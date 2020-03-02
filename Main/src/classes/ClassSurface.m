classdef ClassSurface 
    
   properties            
      bias = {};
      Schottky = 0;
      Ohmic = 0;
      barrier = 0;
   end
    
   methods
      
       function r = CheckValues(Obj)
           r = 0;
           if ~Obj.Ohmic && ~Obj.Schottky
               msgStr = 'error: Could not find the barrier type in ''Surface'' ';
               err    = MException('MATLAB:ValueCheckError',msgStr);
               throw(err);
           end
           r = 1;
       end
   end
   
end