function [Vinitguess,A] = Poisson_Equation_Normal_Compare(OutGenerate,~,~,V0)
% Solves the Normal Poisson Equation
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    0;
Vinit=V0(1);
Vend=V0(end);

% Spatial Distribution of Points
x = OutGenerate.x;

% Number of Nodes for Calculation
No_of_nodes = length(OutGenerate.x);


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

% Starting Value under assumption that n and p are evenly distributed
f=0;
   
 % Construction of b matrix

   % End Boundary Conditions
   % Initial Condition
   b(1,1)=((-1*f*((x(2)-x(1))^2))-Vinit);
   % Final Condition
   b(No_of_nodes-2,1)=((-1*f*((x(No_of_nodes-1)-x(No_of_nodes-2))^2))-Vend);
   % All other Values
   for counter1=2:(No_of_nodes-3)
   b(counter1,1)=(-1*f*((x(counter1+1)-x(counter1))^2));
   end 
   
   
% Solving the Ax = b matrix equation
V=A\b;

% Assigning Values
Vinitguess(1)=V0(1);
Vinitguess(No_of_nodes)=V0(end);

for counter1=1:(No_of_nodes-2)
  Vinitguess(counter1+1)=V(counter1);
end

end