function [H_temp]=getH(OutGenerate)
global Q eps0 

% persistent H H_temp d1 eps1 IsFirst f


x=OutGenerate.x;
epsilon=OutGenerate.epsilon;

No_of_nodes = length(x);



    size_of_matrix = No_of_nodes;
    
    H      = zeros(size_of_matrix+1,size_of_matrix);
    f      = zeros(1,size_of_matrix+1);
    Hdiag  = zeros(1,No_of_nodes-2);
    Hupper = zeros(1,No_of_nodes-3);
    Hlower = zeros(1,No_of_nodes-3);
    d1     = [0 diff(x) 0]; d1(1)   = d1(2); d1(end) = d1(end-1);    
    eps1   = [epsilon(1) epsilon epsilon(end)];
    
    for n=2:No_of_nodes+1
        temp1 = eps1(n-1)/d1(n-1);
        temp2 = eps1(n)/d1(n);
        
        Hdiag(n-1)=(2/(d1(n)+d1(n-1))) *(temp1+temp2);
        
        if n<=No_of_nodes
            Hupper(n-1)=-(2/(d1(n)+d1(n-1)))*temp2;
        end
        
        if n>=3
            Hlower(n-2)=-(2/(d1(n)+d1(n-1)))*temp1;
        end        
    end
    
    H_temp = diag(Hdiag,0)+diag(Hupper,1)+diag(Hlower,-1);

H(1:end-1,:) = H_temp ;
end

