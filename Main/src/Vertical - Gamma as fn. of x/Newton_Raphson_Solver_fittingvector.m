function [ charge ] = Newton_Raphson_Solver_fittingvector( fitting_factor,No_of_nodes,x,charge,Vt,Q,Vbetterguess,ni,n0,p0,indicator )
% Newton Raphson Solver for electron and hole density
  
  % Initialisation of Variables
  % Convergence factor for Inner Newton Raphson 
    alpha       = 1;
  % Maximum Number of Iterations of Inner Newton-Raphson Loop for charge density
    MaxIter     = 600;
  % Tolerence allowed in the charge density profile in the material
    Tolerence   = 1.0E-07;

  % Pre-allocation
    fndash         = zeros(No_of_nodes-2);
    fn             = zeros(No_of_nodes-2,1);
    bn             = zeros(No_of_nodes-2,1);    
    bp             = zeros(No_of_nodes-2,1);
  % Considering only the middle elements of charge
    sqrtdensity       = charge(2:No_of_nodes-1).^0.5;
  
  % Values of Constants
  % b calculation
  bn(end) = fitting_factor(end)/(2*4*1.08*0.5*Q);
  bp(end) = fitting_factor(end)/(2*4*0.57*0.5*Q);
  for i=1:(No_of_nodes-3)
      bn(i) = (fitting_factor(i)+fitting_factor(i+1))/(2*4*1.08*0.5*Q);
      bp(i) = (fitting_factor(i)+fitting_factor(i+1))/(2*4*0.57*0.5*Q);
  end
    phin=(((1/Q)*Vt*log(ni/n0)));
    phip=(((1/Q)*Vt*log(p0/ni)));
     
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
         phi = phip;
        
  elseif indicator == -1
         b   = bn;         
         init= ninit;
         last= nend;
         phi = phin;
  end 
  
  % The Loop solves the Density Gradient Equation for electron and hole density
  iteration = 1;  % Keeps track of the Number of Iterations of the Inner Newton Raphson Loop
   
  while true  
%       fprintf('\n Newton Raphson Loop Iteration: \n');
%       disp(iteration);         

      % Calculation of the electron and hole density equations
      for counter1=1:No_of_nodes-2
          if counter1==1
              fn(counter1,1)=((b(counter1)/((x(1,2)-x(1,1))^2))*(((fitting_factor(counter1+1)/2 +fitting_factor(counter1)/2)*sqrtdensity(counter1+1,1))-(fitting_factor(counter1)+ fitting_factor(counter1+1)/2 )*(sqrtdensity(counter1,1))+(fitting_factor(counter1)/2)*(init)))+(((sqrtdensity(counter1,1))/2)*(Vbetterguess(1,counter1+1)+((indicator)*((1/Q)*Vt*log((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))))-phi));
          elseif counter1==(No_of_nodes-2)
              fn(counter1,1)=((b(counter1)/((x(1,No_of_nodes-1)-x(1,No_of_nodes-2))^2))*((fitting_factor(counter1)*last*0.5)-(fitting_factor(counter1)+fitting_factor(counter1-1)/2 )*(sqrtdensity(counter1,1))+((fitting_factor(counter1)/2+fitting_factor(counter1-1)/2)*sqrtdensity(counter1-1,1))))+(((sqrtdensity(counter1,1))/2)*(Vbetterguess(1,counter1+1)+((indicator)*((1/Q)*Vt*log((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))))-phi));
          else
              fn(counter1,1)=((b(counter1)/((x(1,counter1+1)-x(1,counter1))^2))*((((fitting_factor(counter1+1)/2 +fitting_factor(counter1)/2)*sqrtdensity(counter1+1,1))-(fitting_factor(counter1)+fitting_factor(counter1+1)/2 +fitting_factor(counter1-1)/2)*(sqrtdensity(counter1,1)) +(fitting_factor(counter1)/2 + fitting_factor(counter1-1)/2)*sqrtdensity(counter1-1,1))))+((sqrtdensity(counter1,1)/2)*(Vbetterguess(1,counter1+1)+((indicator)*((1/Q)*Vt*log((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))))-phi));
          end
      end
      % Jacobian Matrix for the electron and hole density equations
      for counter1=1:No_of_nodes-2
          for counter2=1:No_of_nodes-2
              if(counter1==counter2)
                  if counter1 == 1
                  fndash(counter1,counter2)= ((-b(counter1)*(fitting_factor(counter1)+fitting_factor(counter1+1)/2 ))/((x(1,counter1+1)-x(1,counter1))^2)) + ((Vbetterguess(1,counter1+1))/2) - ((phi)/2)  + ((indicator*Vt)/(2*Q))*( 2 + (log(((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))) )  );
                  elseif counter1 == No_of_nodes-2
                  fndash(counter1,counter2)= ((-b(counter1)*(fitting_factor(counter1)+fitting_factor(counter1-1)/2))/((x(1,counter1+1)-x(1,counter1))^2)) + ((Vbetterguess(1,counter1+1))/2) - ((phi)/2)  + ((indicator*Vt)/(2*Q))*( 2 + (log(((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))) )  );    
                  else
                  fndash(counter1,counter2)= ((-b(counter1)*(fitting_factor(counter1)+fitting_factor(counter1+1)/2 +fitting_factor(counter1-1)/2))/((x(1,counter1+1)-x(1,counter1))^2)) + ((Vbetterguess(1,counter1+1))/2) - ((phi)/2)  + ((indicator*Vt)/(2*Q))*( 2 + (log(((sqrtdensity(counter1,1)^2)/ni(1,counter1+1))) )  );   
                  end
              elseif (counter1+1 == counter2)
                  fndash(counter1,counter2) = ((b(counter1)*(fitting_factor(counter1+1)/2 +fitting_factor(counter1)/2))/((x(1,counter1+1)-x(1,counter1))^2));
              elseif (counter1 == counter2+1)
                  fndash(counter1,counter2) = ((b(counter1)*(fitting_factor(counter1)/2 + fitting_factor(counter1-1)/2))/((x(1,counter1+1)-x(1,counter1))^2));
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
      
      norm(sqrtdensitydiff);
      %Printing Iteration number 
      
      % Updating the Variables
      iteration      = iteration+1; 
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