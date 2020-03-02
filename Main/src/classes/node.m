classdef node
    
    properties
        x,y,Na,Nd,epsilon,mn_eff,mp_eff,temperature,
        E_affinity,Eg,Ef,potential,
        layer, sublayer, eta_C, eta_V , Ec, Nc, Nv , Nbt

    end
    
    properties(Dependent = true)
        net_doping,Ev
    end
    
    methods
        
        function net_doping = get.net_doping(Obj)
            net_doping = Obj.Nd - Obj.Na; 
        end 

        function Ev = get.Ev(Obj)
            Ev = Obj.Ec - Obj.Eg; 
        end 
    end
           
    
end