function V_in = Anderson_mixing(V_in,V_out,iter,NoMix,alpha)

global_variables;

persistent VinMat VoutMat derr theta berr IsFirst

if isempty(IsFirst)
    IsFirst = 0;
end

if IsFirst == 0
    VinMat  = zeros(NoMix+1,length(V_in));
    VoutMat = zeros(NoMix+1,length(V_in));
    derr    = zeros(NoMix+1,length(V_in));
    theta   = zeros(NoMix,1);
    berr    = zeros(NoMix,1);
end

IsFirst = IsFirst + 1;

%--------------------------------------------------------------------------

VinMat = [V_in; VinMat ];
VinMat(end,:) = [];

VoutMat = [V_out; VoutMat ];
VoutMat(end,:) = [];

derr = [ V_out - V_in; derr];
derr(end,:) = [];

if iter <= NoMix + 1
    damping_const = 0.001;
    V_in          = V_in*(1-damping_const) + damping_const*(V_out-V_in);
    %     V_in(end) = V_initial(end);
else
    
    Amat = zeros(NoMix,NoMix);
    di = V_out - V_in;
    
    for i=2:NoMix+1    % Note: i must start from 2 only
        for j=2:NoMix+1
            Amat(i-1,j-1) = dot(di-derr(i,:), di-derr(j,:) );
        end
        berr(i-1) = dot(di,di-derr(i,:));
    end
        
    theta = Amat\berr;
    
    VinAvg  = V_in;
    VoutAvg = V_out;
    
    temp = 1 + sum(theta);   % normalizing
    for i=2:NoMix+1    % Note: i must start from 2 only
        VinAvg  = ( VinAvg  + theta(i-1)*(VinMat(i,:) -V_in ) );
        VoutAvg = ( VoutAvg + theta(i-1)*(VoutMat(i,:)-V_out) );
    end
    
    % alpha = 0.05;
    V_in = (1-alpha)*VinAvg + alpha*VoutAvg;
    % V_in      = V_in/abs(V_in(end)) * abs(V_gate);
    % V_in(1)   = V_initial(1);
    
    
end  %end of if

end
