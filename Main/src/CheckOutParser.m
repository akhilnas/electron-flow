function CheckOutParser(OutParser)

if isempty(OutParser.layer)
    msgStr = 'error: No layers found';
    err    = MException('MATLAB:ValueCheckError',msgStr);
    throw(err);
else
    for i =1:length(OutParser.layer)
        OutParser.layer(i).CheckValues();
    end
end

if isempty(OutParser.surface)
    msgStr = 'error: Could not find surface description';
    err    = MException('MATLAB:ValueCheckError',msgStr);
    throw(err);
else
    OutParser.surface.CheckValues();
end

if isempty(OutParser.substrate)
    msgStr = 'error: Could not find substrate description';
    err    = MException('MATLAB:ValueCheckError',msgStr);
    throw(err);
else
    OutParser.substrate.CheckValues();
end

if isempty(OutParser.control)
    msgStr = 'error: Controls field shold not be empty';
    err    =  MException('MATLAB:ValueCheckError',msgStr);
    throw(err);
else 
    OutParser.control.CheckValues(OutParser.layer);
end





end