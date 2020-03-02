%% PARSER_LAYER
% parses layer section of ideal C-V input

%%% function: parser_layer
%  Input: text from input file
% Output: OutParser.layer 
% Called by :
%
% <html> 
% <a href= "parser.html" > "parser.m"  </a> 
% </html>

function [var_layer, NoWords] = parser_layer(C,i,LineNo, Cum_thickness, WORK_DIR)

%% Layer properties:
% <html> 
% <a href= "../../classes/html/ClassLayer.html" > ClassLayer  </a> 
% </html>

var_layer = ClassLayer;
begin     = Cum_thickness;
var_layer.layer_name = C{1}{i+1};

j = 2;

while ((i+j)<= length(C{:}) && ~strcmp(C{1}{i+j},'end'))

    switch C{1}{i+j}

	case 'thickness'
            j = j+1;
            if isnan(str2double(C{1}{i+j}))
                throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
            end
            var_layer.thickness = str2double(C{1}{i+j});
            var_layer.Pos_begin = begin;
			
%       case 'temperature'
%             j = j+1;
%             if isnan(str2double(C{1}{i+j}))
%                 throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
%             end
%             var_layer.temperature = str2double(C{1}{i+j});
                 
% 	  case 'Delta_Ec'
%             j = j+1;
%             if isnan(str2double(C{1}{i+j}))
%                 throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
%             end
%             var_layer.Delta_Ec = str2double(C{1}{i+j});
        case 'grid_spacing'
            j = j+1;
            if isnan(str2double(C{1}{i+j}))
                throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
            end
            var_layer.grid_spacing = str2double(C{1}{i+j});
        case 'Na'
            j = j+1;
            if isnan(str2double(C{1}{i+j}))
                throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
            end
            var_layer.Na = str2double(C{1}{i+j});
        case 'Nd'
            j = j+1;
            if isnan(str2double(C{1}{i+j}))
                throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
            end
            var_layer.Nd = str2double(C{1}{i+j});
        case 'material'
            j = j+1;
            
            if ~exist(fullfile(WORK_DIR,'materials',[C{1}{i+j} '.txt']),'file')
                strMessage = 'error: Unable to find the file ''%s'' (input line No %d)';
                err = MException('MATLAB:FileNotFound',strMessage,[C{1}{i+j} '.txt'],LineNo(i+j));
                throw(err);
            end
            var_layer.material = C{1}{i+j};
			
			%% Material Parsing
		% <html> 
% <a href= "parser_material.html" > Refer parser_material  </a> 
% </html>	
            var_material = parser_material([C{1}{i+j} '.txt'],WORK_DIR);
            var_layer.layer_material = var_material;

            %--------------------------------------------------------------------------------
%% sublayer
%% Layer properties:
% <html> 
% <a href= "../../classes/html/ClassSubLayer.html" > ClassSubLayer  </a> 
% </html>
			case 'sub_layer'
            var_sublayer = ClassSubLayer;
            k = 1;
            while (~strcmp(C{1}{i+j+k},'end'))
                switch C{1}{i+j+k}
                    case 'start'
                        k = k+1;
                        if isnan(str2double(C{1}{i+j+k}))
                            throwErr(C{1}{i+j+k},C{1}{i+j+k-1},LineNo(i+j+k));
                        end
                        var_sublayer.Pos_begin = str2double( C{1}{i+j+k} );
                    case 'thickness'
                        k = k+1;
                        if isnan(str2double(C{1}{i+j+k}))
                            throwErr(C{1}{i+j+k},C{1}{i+j+k-1},LineNo(i+j+k));
                        end
                        var_sublayer.thickness = str2double(C{1}{i+j+k});
                    case 'grid_spacing'
                        k = k+1;
                        if isnan(str2double(C{1}{i+j+k}))
                            throwErr(C{1}{i+j+k},C{1}{i+j+k-1},LineNo(i+j+k));
                        end
                        var_sublayer.grid_spacing = str2double(C{1}{i+j+k});
                    case 'Nd'
                        k = k+1;
                        if isnan(str2double(C{1}{i+j+k}))
                            throwErr(C{1}{i+j+k},C{1}{i+j+k-1},LineNo(i+j+k));
                        end
                        var_sublayer.Nd = str2double(C{1}{i+j+k});
                    case 'Na'
                        k = k+1;
                        if isnan(str2double(C{1}{i+j+k}))
                            throwErr(C{1}{i+j+k},C{1}{i+j+k-1},LineNo(i+j+k));
                        end
                        var_sublayer.Na = str2double(C{1}{i+j+k});
                    case '<#'
                        k = k + comment(C,i+j+k);
                    otherwise
                        msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
                        err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j+k},LineNo(i+j+k));
                        throw(err);
                end
                k = k+1;
            end
            j = j+k;
            var_layer.sublayers = [var_layer.sublayers var_sublayer];
            %--------------------------------------------------------------------------------
        case '<#'
            j = j + comment(C,i+j);
        otherwise
            msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
            err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
            throw(err);
            
    end  %end of switch
    
    j = j+1;
end
NoWords = j;


    function throwErr(value,property,lineNo)
        strMsg1 = 'error: Invalid value ''%s'' detected for property ''%s'' (input line No %d)';
        err1 = MException('MATLAB:InvalidValue',strMsg1,value,property,lineNo);
        throw(err1);
    end
end
