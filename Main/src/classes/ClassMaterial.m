classdef ClassMaterial
    
    properties
        name,
        Eg, 
        mn_eff , mn_t, mn_l,
        mp_eff , mp_lh , mp_hh , mp_so = 0; 
        epsilon, Nc = 0, Nv = 0, E_affinity,
        eta_C = 1, 
        eta_V = 1;  
    end
    
    methods
        
        function mp_eff_val = get.mp_eff(Obj)    
            if      isempty(Obj.mp_eff) && ~isempty(Obj.mp_lh) && ~isempty(Obj.mp_hh)              
                % if mp_eff is not specified, but, mp_lh and mp_hh are
                % specified
                mp_eff_val = ( Obj.mp_lh^1.5 + Obj.mp_hh^1.5 + Obj.mp_so^1.5 )^(2/3);
            elseif  isempty(Obj.mp_eff) && ( isempty(Obj.mp_lh) || isempty(Obj.mp_hh) )             
                % if mp_eff is not specified and, any of mp_lh and mp_hh
                % are not specified
                mp_eff_val = 1;
            else
                % when mp_eff is specified
                mp_eff_val = Obj.mp_eff ;
            end
        end
        
        function mn_eff_val = get.mn_eff(Obj)
            if      isempty(Obj.mn_eff) && ~isempty(Obj.mn_t) && ~isempty(Obj.mn_l)
                % if mn_eff is not specified, but, mn_l and mn_t are
                % specified
                mn_eff_val = (Obj.mn_t * Obj.mn_t * Obj.mn_l)^(1/3);
            elseif  isempty(Obj.mn_eff) && ( isempty(Obj.mn_t) || isempty(Obj.mn_l) ) 
                % if mp_eff is not specified and, any of mn_l and mn_t
                % are not specified
                mn_eff_val = 1;
            else
                % when mn_eff is specified
                mn_eff_val = Obj.mn_eff;
            end
        end
        
        function r = CheckValues(Obj)
            
        end
        
    end
end