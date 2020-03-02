function [ Vbetterguess,A ] = DG_initial( OutGenerate )
%Generates the initial guess for DG Solver

global Q eps0 Kb;
Vt=Kb*OutGenerate.Tempr;

% Initialisation of Variables
Ec               = OutGenerate.Ec;
Ev               = OutGenerate.Ev;
No_of_nodes      = length(OutGenerate.x);

% Pre-allocation of variables/Initial Guess
V0 = zeros(1,No_of_nodes);

BC = struct('B1','','B2','','B1_val',0,'B2_val',0);

if OutGenerate.surface.Schottky || OutGenerate.surface.Ohmic
    BC.B1 = 'Dirichlet';  % BC1 must be dirichlet
end

if OutGenerate.substrate.Schottky || OutGenerate.substrate.Ohmic
    BC.B2 = 'Dirichlet';
elseif OutGenerate.substrate.Zero_Slope
    BC.B2 = 'Neumann';
end

BC.B2='Dirichlet';
% Boundary Condition
V_gate = OutGenerate.surface_potential;
BC.B1_val=V_gate/19.2;
BC.B2_val=V_gate/19.2;

if BC.B1_val==0
    BC.B1_val=1e-13;
end

if strcmp(BC.B1,'Dirichlet')
    V0(1)   = BC.B1_val;
end
if strcmp(BC.B2,'Dirichlet')
    V0(end) = BC.B2_val;
end
disp(V0(end))

        %Normal Poisson Solver
         [Vinitguess , A ]                   = Poisson_Equation_Normal_Compare(OutGenerate,(Ec-Q*V0),(Ev-Q*V0),V0);
        %Self Consistent Poisson Solver
         Vbetterguess                        = Self_Consistent_Poisson_Equation_Compare(OutGenerate,(Ec-Q*Vinitguess),(Ev-Q*Vinitguess),Vinitguess,A,Vt,Q);


end

