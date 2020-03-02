function [ Derivative ] = deln( n,ni,Xold,Vt,bn,phin,Q,deltax )
% Calculates the derivative of n w.r.t potential

   Fn = ( (-phin/(4*n^0.5))  - ( (Vt/(4*Q*n^0.5)) * log(n/ni) ) - (Vt/(2*Q*n^0.5))  + (Xold/(4*n^0.5))  - (bn/((deltax^2)*(n^0.5)))  );
   Derivative  = -((n^0.5)/2)*(1/Fn) ;
end

