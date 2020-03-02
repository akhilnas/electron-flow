%% PARSER
%  Segments input file of ideal c-v
% Calls parsing functions of sections, layer, surface, substrate, controls 

%%% function: parser
%  Input: input file
% Output: Outparser vector
% Called by: 
%
% <html> 
% <a href= "../../html/generate.html" > "generate.m"  </a> 
% </html>


function OutParser = parser(input_file,WORK_DIR)

%% Reading input file

fprintf('Reading the input file...   \n ')

fid = fopen(input_file);

if ~fid
    err = MException('MATLAB:InvalidFileFid',FileName);
    throw(err);
end

C = textscan(fid,'%s','commentStyle', '//');
LineNo = zeros(1,length(C{:}));

frewind(fid);

S = 1; E = 0; CurrentLine = 1;
while ~feof(fid)
    aa = fgetl(fid);
    if ~isempty(aa)
        temp  = textscan(aa,'%s','commentStyle', '//');
        E = S + length(temp{:});
        LineNo(S:E) = CurrentLine;
    end
    CurrentLine = CurrentLine + 1;
    S = E;
end

fclose(fid);

%% Major sections of outparser: layer, surface, substrate, control (Initialising)

OutParser.layer     = [];
OutParser.surface   = [];
OutParser.substrate = [];
OutParser.control  = [];
OutParser.Temperature = 300; % set default temperature to 300. This way even if 'Temeprature' tag is absent from the input file then also temperature gets set to 300K.

i = 1;
CumThickness = 0;
while (i <= length(C{:}))
    switch C{1}{i}
        case '<#'            
            i = i + comment(C,i);
%% Outparser layer 
% <html> 
% <a href= "parser_layer.html" > Refer parser_layer  </a> 
% </html>
			case 'layer'
            [layer,NoWords]     = parser_layer(C,i,LineNo,CumThickness,WORK_DIR);
            i = i + NoWords;            

            OutParser.layer = [OutParser.layer layer];
            CumThickness = CumThickness + layer.thickness;
%%  Outparser surface
% <html> 
% <a href= "parser_surface.html" > Refer parser_surface  </a> 
% </html>
		case 'surface'
		
            [var_surface, NoWords]   = parser_surface(C,i,LineNo);
             i = i + NoWords;
             OutParser.surface = [OutParser.surface var_surface];
        
%% OutParser substrate
% <html> 
% <a href= "parser_substrate.html" > Refer parser_substrate  </a> 
% </html>
		case 'substrate'
            [var_substrate, NoWords] = parser_substrate(C,i,LineNo); 
            i = i + NoWords;
            OutParser.substrate = [OutParser.substrate var_substrate];
%% OutParser controls
% <html> 
% <a href= "parser_controls.html" > Refer parser_controls  </a> 
% </html>
		case 'Controls'
            [var_control, NoWords]  = parser_controls(C,i,LineNo);
            i = i + NoWords;
            OutParser.control = [OutParser.control var_control];
        case 'Temperature'
            if ~isnan(str2double(C{1}{i+1})) && (strcmp(C{1}{i+2},'end'))
                OutParser.Temperature=str2double(C{1}{i+1});
                i=i+2;
            else
                if isnan(str2double(C{1}{i+1}))
                    msgStr = 'error: ''%s'' is not a suitable value as Temperature';
                    err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+1});
                    throw(err);
                end
                if ~(strcmp(C{1}{i+2},'end'))
                    msgStr = 'error: missing "end"';
                    err = MException('MATLAB:InvalidCommand',msgStr);
                    throw(err);
                end
            end
            
        otherwise
            msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
            err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i},LineNo(i));
            throw(err);
    end
    i = i+1;
end

fprintf('Done  \n\n')

end




