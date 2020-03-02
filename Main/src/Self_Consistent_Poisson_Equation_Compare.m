function [Vbetterguess] = Self_Consistent_Poisson_Equation_Compare(OutGenerate,Ec_new,Ev_new,Vinitguess,A,Vt,Q)
  % Solves the Self-Consistent Poisson Equation

  % Initilisation of Variables
  % All the Variable names denote the Convential quantites denoted by these symbols
  epsilon0    = 0.0796.*OutGenerate.epsilon;
  Nd          = OutGenerate.Nd;
  Na          = OutGenerate.Na;
  Ec          = OutGenerate.Ec;
  Ev          = OutGenerate.Ev;
  % Convergence factor for Newton Raphson of Newton-Raphson
  alpha       = 1;
  % Maximum Number of Iterations of Newton-Raphson Loop
  MaxIter     = 5000;
  % Tolerence allowed in the potential profile in the material
  Tolerence   = 1.0E-10;
  
  % Setting Boundary Conditions
  Vinit=Vinitguess(1);
  Vend=Vinitguess(end);

  % Spatial Distribution of Points
  x = OutGenerate.x;

  % Number of Nodes for Calculation
  No_of_nodes = length(OutGenerate.x);
  
  % Pre-allocation
  f=zeros(No_of_nodes-2,1);         % Right-hand side of poisson equation 
  bdash=zeros(No_of_nodes-2);       % Jacobian of AX=b form of the Discrete Poisson equation
  b=zeros(No_of_nodes-2,1);
  Xold=zeros(No_of_nodes-2,1);
  
  % Removing the first and last points out of consideration for the Matrix Solution as they have fixed values
  for counter1=2:(No_of_nodes-1)
      Xold(counter1-1,1) = Vinitguess(counter1);
  end  
  
  % Begining of Newton Raphson Loop for Self-Consistent Poisson Equation
  iteration=1;          % Calculates the Number of Iterations of Newton Raphson Loop

  while iteration>0
    
    % Print the Iteration of the Newton Raphson Loop for the Self-Consistent Poisson Equation    
    %fprintf('Newton Raphson for Self Consistent Poisson Equation = %d \n', iteration);  

    % Charge Density Calculation             
    [~, n_charge, p_charge, ~]=charge_density(Vinitguess,OutGenerate,Ec_new,Ev_new,0,1);
    
    % Calculating the L.H.S of the Poisson Equation
      for counter1=1:(No_of_nodes-2)
        f(counter1,1)=(Q/epsilon0(counter1+1))*((Nd(counter1+1)-Na(counter1+1)+p_charge(counter1+1)-n_charge(counter1+1)));
      end

    % Construction of b matrix

    % End Boundary Conditions
       
       % Initial Condition
       b(1,1)=((-1*f(1,1)*((x(2)-x(1))^2))-Vinit);
       % Final Condition
       b(No_of_nodes-2,1)=((-1*f(No_of_nodes-2,1)*((x(No_of_nodes-1)-x(No_of_nodes-2))^2))-Vend);
       % All other Values
        for counter1=2:(No_of_nodes-3)
          b(counter1,1)=(-1*f(counter1,1)*((x(counter1+1)-x(counter1))^2));
        end
     


    % Construction of Jacobian Matrix
    for counter1=1:(No_of_nodes-2)
        for counter2=1:(No_of_nodes-2)
            if counter1==counter2
                bdash(counter1,counter2)= (-((x(counter1+1)-x(counter1))^2)/epsilon0(counter1+1))*Q*(Q/(Vt))*(((-1)*p_charge(counter1+1))-(n_charge(counter1+1)));
            else
                bdash(counter1,counter2)=0;
            end
        end
    end
    % Newton Raphson formula
    Xdiff =  ((A -bdash)\(A*Xold - b));    
    Xnew  =  Xold - alpha.*Xdiff;
    
    
    
    % Calculating the New Potential Distribution in the material
      for counter1=2:(No_of_nodes-1)
	 	Vinitguess(counter1) = Xnew(counter1-1);
      end
    % Updating Ec and Ev values  
      Ec_new = Ec - Q*Vinitguess;
      Ev_new = Ev - Q*Vinitguess;
      
    % Condition/Tolerence Check 
     if abs(max(Xnew-Xold))<Tolerence
         break;
     end
    % Maximum Number of Iterations of Newton-Raphson Loop Condition Check
    if(iteration>MaxIter)
        break;
    end
    % Updating the Variables
    Xold=Xnew;
    iteration=iteration+1;

  end

  % Assigning the Potential Distribution obtained from the Self-Consistent Poisson Equation
  Vbetterguess(1)=Vinitguess(1);
  Vbetterguess(No_of_nodes)=Vinitguess(end);

    for counter1=1:(No_of_nodes-2)
      Vbetterguess(counter1+1)=Xnew(counter1);
    end

end

