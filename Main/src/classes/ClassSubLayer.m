classdef ClassSubLayer
    
    properties
        node_nos
        thickness
        Pos_begin = 0
        grid_spacing = 0.1
        Nd = 0
        Na = 0
        Nbt = 0
    end
    
    properties(Dependent = true)
        Pos_end
    end
    
    methods
        
        function Pos_end = get.Pos_end(var)
            Pos_end = var.Pos_begin + var.thickness;
        end
        
    end
    
end