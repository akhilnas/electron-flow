%% PARSER_SUBSTRATE
% parses substrate section of ideal C-V input

%%% function: parser_substrate
%  Input: text from input file
% Output: OutParser.substrate 
% Called by :
%
% <html> 
% <a href= "parser.html" > "parser.m"  </a> 
% </html>


function [var_substrate, NoWords] = parser_substrate(C,i,LineNo)

%% ClassSubstrate
% <html> 
% <a href= "../../classes/html/ClassSubstrate.html" > ClassSubstrate  </a> 
% </html>
var_substrate = ClassSubstrate;

j = 1;

while ((i+j)<= length(C{:}) && ~strcmp(C{1}{i+j},'end'))
    switch C{1}{i+j}
        case 'type' 
		%% Type: Schottky, Ohmic, zeroslope
            j = j+1;
            switch C{1}{i+j}
                case 'Schottky'
                    j = j+1;
                    var_substrate.Schottky = 1;
                    if isnan(str2double(C{1}{i+j}))
                       strMsg = 'error: Invalid argument detected for ''Schottky''  (input line No %d)'; 
                       err = MException('MATLAB:InvalidValue',strMsg,LineNo(i+j));
                       throw(err);
                    end
                    var_substrate.barrier = str2double(C{1}{i+j});
                case 'Ohmic'
                    j = j+1;
                    var_substrate.Ohminc = 1;                
                    if isnan(str2double(C{1}{i+j}))
                       strMsg = 'error: Invalid argument detected for ''Ohmic''  (input line No %d)'; 
                       err = MException('MATLAB:InvalidValue',strMsg,LineNo(i+j));
                       throw(err);
                    end
                    var_substrate.barrier = str2double(C{1}{i+j});
                case 'Zero_Slope'                    
                    var_substrate.Zero_Slope = 1;
                otherwise
                    msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
                    err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
                    throw(err);                    
            end
            %% bias
        case 'bias'
             k = 1; A = [];
            while (~strcmp(C{1}{i+j+k},'end') )                
                A = [A ' ' C{1}{i+j+k}];
                k = k+1;
            end
            
            A = textscan(A,'%s','delimiter',';');
            
            for ii = 1:length(A{:})
                var_surface.bias{ii} = [];
                temp = textscan(A{1}{ii},'%s', 'delimiter',',');
                for jj = 1:length(temp{:})
                    temp1 = textscan(temp{1}{jj},'%s');
                    temp2 = str2double(temp1{1});
                    if any(isnan(temp2)) || length(temp2)~=3
                        msgStr = 'error: Invalid arguments found for property ''%s'' (input line No %d) ';
                        err    = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
                        throw(err);
                    end
                    D{ii,jj} = temp2;
                    var_surface.bias{ii} = [var_surface.bias{ii} D{ii,jj}(1):D{ii,jj}(2):D{ii,jj}(3)];
                end                
            end            
            j = j+k;        
        case '<#'
            j = j + comment(C,i+j);            
        otherwise
            msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
            err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j+k},LineNo(i+j+k));
            throw(err);            
    end
    j = j+1;
end    

NoWords = j;    
end

%% TODO: schottky, Ohmic, zero slope? THis version supports only one?