function skip = comment(C,current)
    skip = 1;
    while ~strcmp(C{1}{current+skip},'#>')
        skip = skip+1;    
    end
end
