%% DEFAULT PATH PARSER
% Parses default path in 1DPS_pref.txt

%% 
%%% Function: parser_default_path
%  Input : no argument, but opens 1DPS_pref.txt
% Output: default_path
function default_path = parser_default_path()
 
filename = '1DPS_pref.txt';


fid = fopen(filename);
if fid <= 0
    msgStr = 'error: Could not find the file ''%s'' ';
    err = MException('MATLAB:InvalidFileFid',msgStr,filename);
    throw(err);
end

%% Basic parsing, comment styles etc
C = textscan(fid,'%s','CommentStyle','#','delimiter','= ','MultipleDelimsAsOne','1');
fclose(fid);

C = strtrim(C);

i = 1;
%% Checks for 'default_path' string 
while ~strcmp(C{1}{i},'default_path')
    if i == length(C{1})
        msgStr = 'error: Invalid command found in ''%s''';
        err = MException('MATLAB:InvalidCommand',msgStr,filename);
        throw(err);
    end
    i = i+1;
end
i = i+1;
%% checks if text exists, and assigns the output 'default_path' to it 
if (length(C{1}) -i + 1) > 1
    msgStr = 'error: Invalid command found in ''%s''';
    err = MException('MATLAB:InvalidCommand',msgStr,filename);
    throw(err);    

     
else
    temp = textscan(C{1}{i},'%q');
    default_path = temp{1}{1};
end

if ~exist(default_path,'dir')
    msgStr = 'error: Could not find the folder ''%s'' ';
    err = MException('MATLAB:FolderNotExist',msgStr,default_path);
    throw(err);
end

end