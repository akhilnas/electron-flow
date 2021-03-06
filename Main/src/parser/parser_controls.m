%% PARSER_CONTROLS
%  parses controls section of ideal C-V input
% Output directory parsing for ideal C-V (attaches file prefix)
% Border trap solvers

%%% function: parser_controls
%  Input: text from input file
% Output: OutParser.controls 
% Called by :
%
% <html> 
% <a href= "parser.html" > "parser.m"  </a> 
% </html>





function [var_Controls, NoWords] = parser_controls(C,i,LineNo)

var_Controls = ClassControls;
%% ClassControls
% <html> 
% <a href= "../../classes/html/ClassControls.html" > ClassControls  </a> 
% </html>

j = 1;

while ((i+j)<= length(C{:}) && ~strcmp(C{1}{i+j},'end'))
    switch C{1}{i+j}
	    %% Solvers
        case 'Solvers'
            j = j+1;
            while (~strcmp(C{1}{i+j},'end'))
                switch C{1}{i+j}
                 %* Schrodinger POisson   
				  case 'Schrodinger_Poisson'
                        j = j+1;
                        var_Controls.Solvers.SP = 1;
                        while (~strcmp(C{1}{i+j},'end'))
                            switch C{1}{i+j}
                       			case 'SchrodingerStart'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.SchrodingerStart = str2double(C{1}{i+j});
                                case 'SchrodingerStop'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.SchrodingerStop  = str2double(C{1}{i+j});
                                case 'max_iterations'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.max_iterations   = str2double(C{1}{i+j});
                                case 'tolerance'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.tolerance  = str2double(C{1}{i+j});
                                %% Mixing parameters
								case 'mixing'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.mixing   = str2double(C{1}{i+j});
                                case 'alpha'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.alpha    = str2double(C{1}{i+j});
                                case '<#'
                                    j = j + comment(C,i+j);
                                otherwise
                                    msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
                                    err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
                                    throw(err);
                            end
                            j = j+1;
                        end
                    j=j+1;    
                    case 'Poisson'
                        j = j+1;
                        var_Controls.Solvers.P = 1;
                        while (~strcmp(C{1}{i+j},'end'))
                            switch C{1}{i+j}
                       			
                                case 'max_iterations'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.max_iterations   = str2double(C{1}{i+j});
                                case 'tolerance'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.tolerance  = str2double(C{1}{i+j});
                                %% Mixing parameters
								case 'mixing'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.mixing   = str2double(C{1}{i+j});
                                case 'alpha'
                                    j = j+1;
                                    if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                    end
                                    var_Controls.Schodinger_Poisson.alpha    = str2double(C{1}{i+j});
                                case '<#'
                                    j = j + comment(C,i+j);
                                otherwise
                                    msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
                                    err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
                                    throw(err);
                            end
                            j = j+1;
                        end
                        j = j+1;
						%% Border Trap solving
                    case 'Border_Trap'
                        j = j+1;
                        var_Controls.Solvers.Trap = 1;
                        while (~strcmp(C{1}{i+j},'end'))
                           switch C{1}{i+j}
                               case 'Trap_Layer'
                                   var_Controls.Border_Trap.Trap_Layer = C{1}{i+j};  
                               case 'Injection_Layer'
                                   var_Controls.Border_Trap.Injecting_Layer = C{1}{i+j}; 
                               case 'grid_spacing'
                                   if isnan(str2double(C{1}{i+j}))
                                        throwErr(C{1}{i+j},C{1}{i+j-1},LineNo(i+j));
                                   end
                                   var_Controls.Border_Trap.grid_spacing = str2double(C{1}{i+j});
                               case '<#'
                                   j = j + comment(C,i+j);
                               otherwise
                                   msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
                                   err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
                                   throw(err);
                           end
                           j = j+1;
                        end                        
                    otherwise
                        msgStr = 'error: Unable to evaulate the command ''%s'' (input line No %d)';
                        err = MException('MATLAB:InvalidCommand',msgStr,C{1}{i+j},LineNo(i+j));
                        throw(err);
                end
            end
			%% Output directory code
        case 'Output_directory'
            temp1 = j;
            j = j+1;
            temp = {};
            while (~strcmp(C{1}{i+j},'end'))
                if strcmp(C{1}{i+j},'<#')                    
                    j = j + comment(C,i+j);
                else
                    temp = {temp{:} C{1}{i+j}};
                end
                j = j+1;
            end
            if length(temp) ~= 1
                msgStr = 'error: Invalid arguments found for ''Output_directory'' (input line No %d)';
                err = MException('MATLAB:InvalidCommand',msgStr,LineNo(i+temp1));
                throw(err);
            end
            var_Controls.Output_directory = temp{1} ;
            % Create the directory to make sure the directory can be
            % created with the given path and name
            if ~exist(var_Controls.Output_directory,'dir')
                [status,message] = mkdir(var_Controls.Output_directory);
                if ~status
                    msgStr = 'error: Can not create the directory ''%s'' (input line No %d): %s ';
                    err = MException('MATLAB:MKDIR:OSError',msgStr,var_Controls.Output_directory,LineNo(temp1),message);
                    throw(err);
                end
            else 
                fid = fopen(fullfile(var_Controls.Output_directory,'test'),'w');
                if fid <=0
                    msgStr = 'error: Folder ''%s'' is write protected (input line No %d) ';
                    err = MException('MATLAB:MKDIR:OSError',msgStr,var_Controls.Output_directory,LineNo(i+temp1));
                    throw(err);
                else 
                    delete(fullfile(var_Controls.Output_directory,'test'));
                end
            end
          %% file_prefix  
        case 'file_prefix'
            temp1 = j;
            j = j+1;
            temp = {};
            while (~strcmp(C{1}{i+j},'end'))
                if strcmp(C{1}{i+j},'<#')
                    j = j + comment(C,i+j);
                else
                    temp = {temp{:} C{1}{i+j}};
                end
                j = j+1;
            end
            if length(temp) ~= 1
                msgStr = 'error: Invalid arguments found for ''file_prefix'' (input line No %d)';
                err = MException('MATLAB:InvalidCommand',msgStr,LineNo(i+temp1));
                throw(err);
            end
            var_Controls.file_prefix = temp{:} ;
            clear temp temp1;
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
%% Error for controls and output dir
    function throwErr(value,property,lineNo)
        strMsg1 = 'error: Invalid value ''%s'' detected for property ''%s'' (input line No %d)';
        err1 = MException('MATLAB:InvalidValue',strMsg1,value,property,lineNo);
        throw(err1);
    end

end
%% TODO: Border trap dependence?