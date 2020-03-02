function [ V ] = Poisson( n,p,OutGenerate )
%Solves Poisson Equation

global Q eps0 Kb ;


% Assigning
epsilon0  = 0.0796.*OutGenerate.epsilon;
Nd        = OutGenerate.Nd;
Na        = OutGenerate.Na;
No_of_nodes = length(OutGenerate.x);
Vinit    = OutGenerate.surface_potential/19.2;
Vend     = OutGenerate.surface_potential/19.2;
V        = zeros(1,No_of_nodes);
x        = OutGenerate.x;
% Calculating the L.H.S of the Poisson Equation
      for counter1=1:(No_of_nodes-2)
        f(counter1,1)=(Q/epsilon0(counter1+1))*((Nd(counter1+1)-Na(counter1+1)+p(counter1+1)-n(counter1+1)));
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
        
    % Construction of the A matrix

   % Implementing first Boundary Condition
   A(1,1)=-2;
   A(1,2)=1;
   for count=3:(No_of_nodes-2)
       A(1,count)=0;
   end
   
   
 % Implementing the rest of the Matrix
 
 for counter1=2:(No_of_nodes-3)
     for counter2=1:(No_of_nodes-2)
         if counter2==(counter1-1)
         A(counter1,counter2)=1;
         elseif counter2==counter1
         A(counter1,counter2)=-2;
         elseif counter2==(counter1+1)
         A(counter1,counter2)=1;
         else
         A(counter1,counter2)=0;
         end
     end
 end
 
  % Implementing Last Boundary Condition
   A(No_of_nodes-2,No_of_nodes-2)=-2;
   A(No_of_nodes-2,No_of_nodes-3)=1;
   for count=1:(No_of_nodes-4)
       A(No_of_nodes-2,count)=0;
   end
   
  % Saving on Memory
   A = sparse(A); 
   
   %Solving
   
   x = A\b;
   
   %Profile 
   V(1) = Vinit;
   V(end) = Vend;
   V(2:No_of_nodes-1) = x;
     



end

