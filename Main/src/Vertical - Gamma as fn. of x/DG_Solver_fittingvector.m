function [VDG , n_charge_DG , p_charge_DG] = DG_Solver_fittingvector(fitting_factor,OutGenerate,Vbetterguess,A,Vt,Q)
  % Main part of Density Gradient Implementation
  % Solver to find the Consistent solution to the poisson and density gradient equations for hole and electron density

  % Initialisation of Variables
  % All the Variable names denote the Convential quantites denoted by these symbols
    epsilon0    = 0.0796.*OutGenerate.epsilon;
    Nd          = OutGenerate.Nd;
    x           = OutGenerate.x;
    No_of_nodes = length(OutGenerate.x);
    Nc          = OutGenerate.Nc;
    Nv          = OutGenerate.Nv;
    Ec          = OutGenerate.Ec;
    Ev          = OutGenerate.Ev;
    Ef          = OutGenerate.Ef;

  % Calculating intrinsic charge carrier concentration and n0 values
    ni          = sqrt(Nc.*Nv.*exp((1/Vt)*(Ev-Ec)));
    n0          = Nc.*exp((1/Vt)*(Ef-Ec));
    p0          = (ni.^2)./n0;
  
  % Convergence factor for Newton Raphson of DG Newton-Raphson
    alpha       =   1;
  % Maximum Number of Iterations of DG Newton-Raphson Loop for potential
    MaxIter     = 500;
  % Tolerence allowed in the potential profile in the material
    Tolerence   = 1.0E-15;
  % Pre-allocation  
    ffinal      = zeros(No_of_nodes-2,1);
    bdash       = zeros(No_of_nodes-2);
  % Construction of Potential function under consideration
    Xold        = Vbetterguess(2:No_of_nodes-1)';
    
    
  % Initial Charge Density Calculation
    % Charge Density Calculation             
        [~, n_charge_DG, p_charge_DG, ~] = charge_density(Vbetterguess,OutGenerate,(Ec-Q*Vbetterguess),(Ev-Q*Vbetterguess),1);
    
  % Begining of Outer Newton Raphson Loop for the Poisson Equation
  % The Outer Loop solves the Poisson Equation
  % The Inner Loops solve the Density Gradient Equation for electron and hole density
	iteration=1;  % Keeps track of the Number of Iterations of the Outer Newton Raphson Loop   
    

    while true    
	  
% 	    disp(iteration);
	    
	    % Charge Density Calculation             
%         [~, n_charge_DG, p_charge_DG, ~] = charge_density(Vbetterguess,OutGenerate,(Ec-Q*Vbetterguess),(Ev-Q*Vbetterguess),1);
	      
        % Taking the Transpose for ease of Calculation in the following functions
        n_charge_DG = n_charge_DG' ;
        p_charge_DG = p_charge_DG' ; 
%         
%         n_charge_DG = zeros(No_of_nodes-2,1);
         	  
	    % Solving the Density Gradient Equation for electron density
	    [ n_charge_DG ] = Newton_Raphson_Solver_fittingvector( fitting_factor,No_of_nodes,x,n_charge_DG,Vt,Q,Vbetterguess,ni,n0,p0,-1 );

	    % Solving the Density Gradient Equation for hole density
	    [ p_charge_DG ] = Newton_Raphson_Solver_fittingvector( fitting_factor,No_of_nodes,x,p_charge_DG,Vt,Q,Vbetterguess,ni,n0,p0,1 );

	    % Construcion of b matrix of poisson equation
        for counter1=1:No_of_nodes-2
	        if counter1==1
	             ffinal(counter1,1)=(((x(2)-x(1))^2)/epsilon0(counter1+1))*Q*((p_charge_DG(counter1+1))-(n_charge_DG(counter1+1))+Nd(counter1+1))+ Vbetterguess(1);
	        elseif counter1==(No_of_nodes-2)
	             ffinal(counter1,1)=(((x(No_of_nodes-1)-x(No_of_nodes-2))^2)/epsilon0(counter1+1))*Q*((p_charge_DG(counter1+1))-(n_charge_DG(counter1+1))+Nd(counter1+1)) + Vbetterguess(No_of_nodes);
	        else
	             ffinal(counter1,1)=(((x(counter1+1)-x(counter1))^2)/epsilon0(counter1+1))*Q*((p_charge_DG(counter1+1))-(n_charge_DG(counter1+1))+Nd(counter1+1));
	        end
        end
        
       % Construction of Jacobian Matrix of the Poisson Equation (Boltzmann Approximation)
            for counter1=1:(No_of_nodes-2)
               for counter2=1:(No_of_nodes-2)
                if counter1==counter2
                    bdash(counter1,counter2)= (((x(counter1+1)-x(counter1))^2)/(epsilon0(counter1+1)))*Q*(Q/Vt)*((-p_charge_DG(counter1+1))-(n_charge_DG(counter1+1)));
                else
                    bdash(counter1,counter2)=0;
                end
               end
            end
        
        % Newton Raphson formula        
	    Xdiff=((A + bdash)\(A*Xold + ffinal));
        Xold=Xold-(alpha.*Xdiff);
    
        % Calculating the New Potential Distribution in the material
        for counter1=2:(No_of_nodes-1)
	 	Vbetterguess(counter1) = Xold(counter1-1);
        end
        
        % Condition/Tolerence Check 	  
        if max(abs(Xdiff))<Tolerence
	     break;
        end
        
        % Checking for NaN Error
        if isnan(Xdiff)
           break;
        end
        norm(Xdiff);
        % Maximum Number of Iterations of Outer Newton-Raphson Loop Condition Check
        if(iteration>MaxIter)
	     break;
        end

	    % Updating the Variables
        iteration=iteration+1;
    end
	

	% Construction of Potential Profile obtained from Density Gradient Model
	VDG(1)           = Vbetterguess(1);
	VDG(No_of_nodes) = Vbetterguess(No_of_nodes);
    for counter1=2:(No_of_nodes-1)
            VDG(counter1)= Xold(counter1-1);
    end

end