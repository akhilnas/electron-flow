function [ charge ] = Newton_Raphson_Solver1( fitting_factor,No_of_nodes,x,charge,Vt,Q,Vbetterguess,ni,Nc, indicator )
% Newton Raphson Solver for electron and hole density
  
  % Initialisation of Variables
  % Convergence factor for Inner Newton Raphson 
    alpha       = 1;
  % Maximum Number of Iterations of Inner Newton-Raphson Loop for charge density
    MaxIter     = 60000;
  % Tolerence allowed in the charge density profile in the material
    Tolerence   = 5.0E-05;

  % Pre-allocation
    fndash         = zeros(No_of_nodes-2);
    fn             = zeros(No_of_nodes-2,1);  
    
  % Considering only the middle elements of charge
    sqrtdensity       = charge(2:No_of_nodes-1).^0.5;
    N_middle          = Nc(2:No_of_nodes-1)';
  
  % Values of Constants
    bn  = fitting_factor/(4*sqrt(2)*0.5*1.08);            %bn calculation
    bp  = fitting_factor/(12*0.57*Vt);            %bp calculation
    phi = Vt*log((sqrtdensity.^2)./N_middle);
     
  % Electron and hole density boundary conditions
    ninit = 0;
    nend  = 0;
    pinit = 0;    
    pend  = 0;

  % Selecting the Constants depending on the sign of chrage Density that we want to solve for  
  if indicator == 1
         b   = bp;         
         init= pinit;
         last= pend;
        
  elseif indicator == -1
         b   = bn;         
         init= ninit;
         last= nend;
  end 
  
  % The Loop solves the Density Gradient Equation for electron and hole density
  iteration = 1;  % Keeps track of the Number of Iterations of the Inner Newton Raphson Loop
   
  while true  
      fprintf('\n Newton Raphson Loop Iteration: \n');
      disp(iteration);         

      % Calculation of the electron and hole density equations
      for counter1=1:No_of_nodes-2
          if counter1==1
              fn(counter1,1)=((b/((x(1,2)-x(1,1))^2))*((sqrtdensity(counter1+1,1))-2*(sqrtdensity(counter1,1))+(init)))+(((sqrtdensity(counter1,1))/2)*(Vbetterguess(1,counter1+1)+((indicator)*((1/Q)*Vt*log((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))))-phi(counter1,1)));
          elseif counter1==(No_of_nodes-2)
              fn(counter1,1)=((b/((x(1,No_of_nodes-1)-x(1,No_of_nodes-2))^2))*((last)-2*(sqrtdensity(counter1,1))+(sqrtdensity(counter1-1,1))))+(((sqrtdensity(counter1,1))/2)*(Vbetterguess(1,counter1+1)+((indicator)*((1/Q)*Vt*log((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))))-phi(counter1,1)));
          else
              fn(counter1,1)=((b/((x(1,counter1+1)-x(1,counter1))^2))*((sqrtdensity(counter1+1,1))-2*(sqrtdensity(counter1,1)) +(sqrtdensity(counter1-1,1))))+((sqrtdensity(counter1,1)/2)*(Vbetterguess(1,counter1+1)+((indicator)*((1/Q)*Vt*log((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))))-phi(counter1,1)));
          end
      end
      % Jacobian Matrix for the electron and hole density equations
      for counter1=1:No_of_nodes-2
          for counter2=1:No_of_nodes-2
              if(counter1==counter2)
                  fndash(counter1,counter2)= ((-b*2)/((x(1,counter1+1)-x(1,counter1))^2)) + ((Vbetterguess(1,counter1+1))/2) - ((phi(counter1,1))/2) - Vt + ((indicator*Vt)/(2*Q))*( 2 + (log(((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))) )  );                          
              elseif (counter1+1 == counter2)
                  fndash(counter1,counter2) = ((b)/((x(1,counter1+1)-x(1,counter1))^2));
              elseif (counter1 == counter2+1)
                  fndash(counter1,counter2) = ((b)/((x(1,counter1+1)-x(1,counter1))^2));
              else
                  fndash(counter1,counter2) = 0;
              end
          end
      end
      fndash = sparse(fndash);
      % Newton Raphson formula
      sqrtdensitydiff = (fndash\fn);
      sqrtdensity     = sqrtdensity - (alpha.*sqrtdensitydiff);
              
      % Condition/Tolerence Check      
      if max(abs(sqrtdensitydiff))<Tolerence
          break;
      end
      
      % Checking for NaN Error
      if isnan(sqrtdensitydiff)
          break;
      end
      
      % Maximum Number of Iterations of Inner Newton-Raphson Loop Condition Check
      if(iteration>MaxIter)
          break;
      end     
      
      norm(sqrtdensitydiff)
      %Printing Iteration number 
      
      % Updating the Variables
      iteration      = iteration+1; 
      phi            = Vt*log((sqrtdensity.^2)./N_middle);
  end

  % Reassigning
  charge(1)           = init;
  charge(No_of_nodes) = last;
  for counter1=1:(No_of_nodes-2)
    charge(counter1+1) = sqrtdensity(counter1)^2;
  end

  % Reorienting the Matrix
  charge = charge' ;

end

