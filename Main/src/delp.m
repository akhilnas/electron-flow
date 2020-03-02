function [ Derivative ] = delp( p,ni,Xold,Vt,bp,phip,Q,deltax )
% Calculates the derivative of n w.r.t potential

   Fn = ( (-phip/(4*p^0.5))  + ( (Vt/(4*Q*p^0.5)) * log(p/ni) ) + (Vt/(2*Q*p^0.5))  + (Xold/(4*p^0.5))  - (bp/((deltax^2)*(p^0.5)))  );
   Derivative  = -((p^0.5)/2)*(1/Fn) ;
end

