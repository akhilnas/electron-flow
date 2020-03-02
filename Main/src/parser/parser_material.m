%% PARSER_MATERIAL
% parses material properties

%%% function: parser_material
%  Input: material file
% Output: ClassMaterial; used in OutParser.layer 
% Called by:
%
% <html> 
% <a href= "parser_layer.html" >  "parser_layer" </a> 
% </html>

function var_material = parser_material(filename, WORK_DIR )

     % Materials folder in default_path (stored in WORK_DIR)
     pathname = [ WORK_DIR '/materials/'];
         
     fid = fopen([pathname filename]);
       D = textscan(fid,'%s','commentStyle', '//');
     fclose(fid);
%% ClassMaterial
% <html> 
% <a href= "../../classes/html/ClassMaterial.html" > ClassMaterial  </a> 
% </html>
%% TODO: which properties are compulsory?
     var_material = ClassMaterial;
     p = 1;     
     while (p <= length(D{:}))         
         switch D{1}{p}
		 % Name
             case 'name'
                 p = p+1;
                 var_material.name = D{1}{p};
             
			 case 'Eg'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                    throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.Eg = str2double(D{1}{p});
             
			 case 'mn'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.mn_eff = str2double(D{1}{p});
             
			 case 'mn_t'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.mn_t = str2double(D{1}{p});                 
             
			 case 'mn_l'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.mn_l = str2double(D{1}{p});
             
			 case 'mp'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.mp_eff = str2double(D{1}{p});
             
			 case 'mp_lh'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.mp_lh = str2double(D{1}{p});                 
             
			 case 'mp_hh'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.mp_hh = str2double(D{1}{p});
             
			 case 'epsilon'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.epsilon = str2double(D{1}{p});
             
			 case 'Nc'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.Nc = str2double(D{1}{p});
             
			 case 'Nv'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.Nv = str2double(D{1}{p});
             
			 case 'E_affinity'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.E_affinity = str2double(D{1}{p});
             
			 case 'eta_C'
                 p = p+1;
                 if (isnan(str2double(D{1}{p})))
                     throwErr(D{1}{p},D{1}{p-1},filename);
                 end
                 var_material.eta_C = str2double(D{1}{p});

             otherwise
                  msgStr = 'error: Unable to evaulate the command ''%s'' (in file ''%s'')';
                  err = MException('MATLAB:InvalidCommand',msgStr,D{1}{p},filename);
                  throw(err);
         end
         p = p+1;
     end

	 % Error format
    function throwErr(value, property, filename)
        msgStr = 'error: Invalid value ''%s'' detected for property ''%s'' (in file ''%s'')';
        err = MException('MATLAB:InvalidValue',msgStr,value,property,filename);
        throw(err);
    end
     
end