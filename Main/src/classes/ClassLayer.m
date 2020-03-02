classdef ClassLayer
    
    properties
        layer_name
        nodes
        thickness
        Pos_begin = 0
        material
        Delta_Ec = 0
        Ec
        temperature = 300
        sublayers
        grid_spacing = 0.1
        Nd = 0
        Na = 0
        Nbt = 0
        
        layer_material = ClassMaterial;
    end
    
    properties(Dependent = true)
        Pos_end
    end
    
    methods
        
        function Pos_end = get.Pos_end(Obj)
            Pos_end = Obj.Pos_begin + Obj.thickness;
        end
        
        function r = CheckValues(Obj)
            r = 0;
            t_layer    = Obj.thickness;
            t_sublayer = 0;
            if ~isempty(Obj.sublayers)
                for j = 1:length(Obj.sublayers)
                    
                    if Obj.sublayers(j).Pos_begin > t_layer
                        msgStr = 'error: ''start'' must be < layer ''thickness'' (in layer ''%s'', sublayer ''%d'')';
                        err    = MException('MATLAB:ValueCheckError',msgStr,Obj.layer_name,j);
                        throw(err);
                    end
                    t_sublayer = t_sublayer + Obj.sublayers(j).thickness;
                end
                if t_layer < t_sublayer
                    msgStr = 'error: Sum of the thicknesses of sublayers > thickness of the layer (in layer ''%s'')';
                    err    = MException('MATLAB:ValueCheckError',msgStr,Obj.layer_name);
                    throw(err);
                end
            end
            r = 1;
        end
                
        
    end
    
end