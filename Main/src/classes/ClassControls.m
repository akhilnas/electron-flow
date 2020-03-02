%% ClassControls
% 

classdef ClassControls
    
    properties
        Solvers = struct('SP',0,'Trap',0,'IV',0)
        Output_directory = '/home/Output/'
        Schodinger_Poisson = struct('SchrodingerStart',[],'SchrodingerStop',[], ...
            'max_iterations',100,'tolerance',1e-4,'mixing',2,'alpha',0.01)
        Border_Trap = struct('Trap_Layer',[],'Trap_Layer_no',[],'Injecting_Layer',[],'Injecting_Layer_no',[],'grid_spacing',[])
        file_prefix = 'out_'
    end
    
    methods
        
        function r = CheckValues(Obj,layer)
            r = 0;
            if Obj.Solvers.SP
                if  isempty(Obj.Schodinger_Poisson.SchrodingerStart)
                    msgStr = 'error: Could not find value of property ''SchrodingerStart''';
                    err = MException('MATLAB:ValueNotSet',msgStr);
                    throw(err);
                elseif isempty(Obj.Schodinger_Poisson.SchrodingerStop)
                    msgStr = 'error: Could not find value of property ''SchrodingerStop''';
                    err = MException('MATLAB:ValueNotSet',msgStr);
                    throw(err);
                end
                
                if ~Obj.Solvers.SP && ~Obj.Solvers.IV
                    msgStr = 'error: No solvers specified in ''Controls'' field';
                    err = MException('MATLAB:ValueNotSet',msgStr);
                    throw(err);
                end
            end
            if Obj.Solvers.Trap
                temp = false;
                for i=1:length(layer) 
                    if strcmp(layer(i).layer_name,'Trap_Layer')
                        temp = true;
                    end
                end
                if ~temp
                    msgStr = 'error: Invalid value found for property ''Trap_Layer'' ';
                    err = MException('MATLAB:InvalidValue',msgStr);
                    throw(err);
                end
                temp = false;
                for i=1:length(layer) 
                    if strcmp(layer(i).layer_name,'Injection_Layer')
                        temp = true;
                    end
                end
                if ~temp
                    msgStr = 'error: Invalid value found for property ''Injection_Layer'' ';
                    err = MException('MATLAB:InvalidValue',msgStr);
                    throw(err);
                end
            end
            r = 1;
        end
        
        function layer_no = layer_no(Obj,layer)
            if isempty(Obj.Border_Trap.Trap_Layer_no)
                for i=1:length(layer) 
                    if strcmp(layer(i),Obj.Border_Trap.Trap_Layer)
                        Obj.Border_Trap.Trap_Layer_no = i;
                        layer_no = Obj.Border_Trap.Trap_Layer_no;
                    end
                end
            else
                layer_no = Obj.Border_Trap.Trap_Layer_no;
            end
        end
        
        function layer_no = Injecting_Layer_no(Obj,layer)
            if isempty(Obj.Border_Trap.Injecting_Layer_no)
                for i=1:length(layer) 
                    if strcmp(layer(i),Obj.Border_Trap.Injecting_Layer)
                        Obj.Border_Trap.Injecting_Layer_no = i;
                        layer_no = Obj.Border_Trap.Injecting_Layer_no;
                    end
                end
            else
                layer_no = Obj.Border_Trap.Injecting_Layer_no;
            end
        end
                
    end
    
end