%% PARSER_SURFACE
% parses surface section of ideal C-V input

%%% function: parser_surface
%  Input: text from input file
% Output: OutParser.surface 
% Called by :
%
% <html> 
% <a href= "parser.html" > "parser.m"  </a> 
% </html>

function [var_surface, NoWords] = parser_surface(C,i,LineNo)

%% ClassSurface
% <html> 
% <a href= "../../classes/html/ClassSurface.html" > ClassSurface  </a> 
% </html>
var_surface = ClassSurface;

j = 1;

while ((i+j)<= length(C{:}) && ~strcmp(C{1}{i+j},'end'))
    switch C{1}{i+j}
	    %% Type: Schottky, Ohmic
        case 'type'
            j = j+1;
            switch C{1}{i+j}
                case 'Schottky'
                    j = j+1;
                    var_surface.Schottky = 1;
                    if isnan(str2double(C{1}{i+j}))
                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                    end
                    var_surface.barrier = str2double(C{1}{i+j});
                    fprintf('var_surface.barrier')
                case 'Ohmic'
                    j = j+1;
                    var_surface.Ohminc = 1;
                    if isnan(str2double(C{1}{i+j}))
                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                    end
                    var_surface.barrier = str2double(C{1}{i+j});
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
            err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
            throw(err);
    end
    j = j+1;
end

NoWords = j;

    function throwErr(value,property,lineNo)
        strMsg1 = 'error: Invalid value ''%s'' detected for property ''%s'' (input line No %d)';
        err1 = MException('MATLAB:InvalidValue',strMsg1,value,property,lineNo);
        throw(err1);
    end

end
%% TODO: schottky and ohmic surface?