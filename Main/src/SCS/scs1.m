% This file is Newton Raphson implementation of CVsimulator.
% HFCV versions
% Status - Developed
% Comments - With itegration of this code poisson.m will become deprecated.
% New function getH() will be added to src. Potential is no more treated as
% energy. E=qV, where q=sqrt(2) in Rydberg units. A new control keyword
% 'Method' NR/NR+Anderson will get introduced. For now drho/dV is less
% accurate in SP. Need to implement more accurate form.
%
%
% Created by : Dhirendra (dhirendra22121987@gmail.com)
% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

function scs(H,OutGenerate,V_ramp)


global Q eps0 Kb;
Vt=Kb*OutGenerate.Tempr;

'here'

%OutGenerate=getEquilibrium(H,OutGenerate);


max_iter         = OutGenerate.control.Schodinger_Poisson.max_iterations;
tolerance        = OutGenerate.control.Schodinger_Poisson.tolerance;
output_dir       = OutGenerate.control.Output_directory;
file_prefix      = OutGenerate.control.file_prefix;
NoMix            = OutGenerate.control.Schodinger_Poisson.mixing;
alpha            = OutGenerate.control.Schodinger_Poisson.alpha;
Phi_ms           = OutGenerate.Phi_ms;
Ec               = OutGenerate.Ec;
Ev               = OutGenerate.Ev;
No_of_nodes      = length(OutGenerate.x);
x                = OutGenerate.x;
epsilon          = OutGenerate.epsilon;
N=length(x);

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
% Boundary conditions on H
% Symmetry lost in this operations. Need to find a way to make H symmetric.
H(1,:)=0;
H(1,1)=(epsilon(1)/((x(2)-x(1))^2));
% Hend=H(end,:);
H(end,:)=0;
H(end,end)=(epsilon(end)/((x(end)-x(end-1))^2));
%H(end,end-1)=-(epsilon(end)/((x(end)-x(end-1))^2))/2;
%H(end,end-1)=-(epsilon(end)/((x(end)-x(end-1))^2))/2;



% V_in  = zeros(1,No_of_nodes);
% V_out = zeros(1,No_of_nodes);
V0 = zeros(1,No_of_nodes);
% Phi_ms=0.5; % for now Phi_ms is not used


out_dir{1}=OutGenerate.control.Output_directory;
if output_dir(end)=='/'
    if ~exist([output_dir(1:end-1) '_HF'],'dir')
        mkdir([output_dir(1:end-1) '_HF']);
    end
    out_dir{2}=[output_dir(1:end-1) '_HF'];
else
    mkdir([output_dir '_HF']);
    out_dir{2}=[output_dir '_HF'];
end


V_gate = 0.5;
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
%         if aa==1
%             NoMix            = OutGenerate.control.Schodinger_Poisson.mixing;
%             alpha            = OutGenerate.control.Schodinger_Poisson.alpha;
%         end
%         if aa==2
%             NoMix            = 0
%             alpha            = 1
%         end

iter=0;
check=1;

disp(OutGenerate.x(end))
while(check)
    iter=iter+1;
    
    if(OutGenerate.control.Solvers.SP)
        OutSchr = Schrodinger4(OutGenerate,(Ec-Q*V0),(Ev-Q*V0));
        OutSchr1 = Schrodinger4(OutGenerate,(Ec+0.0132-Q*V0),(Ev-Q*V0));
        [rho, n_charge, p_charge, rho1] = charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0),OutSchr);
        
    else
        [rho n_charge p_charge rho1]=charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0));
    end
    
    rho=Q*(OutGenerate.Nd-OutGenerate.Na+p_charge-n_charge);
    
    dn=(Q/Vt)*n_charge;
    dp=(-Q/Vt)*p_charge;
    
    f=rho'/eps0;
    drho=Q*(dp-dn);
    df=drho/eps0;
    NR=H;
    
    for j=2:N-1
        NR(j,j)=NR(j,j)-df(j);
    end
    
    
    
    R=H*V0'-f;
    R(1)=0;
    R(end)=0;
    
    dV=-NR\R;
    
    V0=V0+dV';

    %V0(1)=BC.B1_val;
    %V0(end)=BC.B2_val;
    
    fprintf('Iteration = %d \t Bias = %.2f \t Error in V = %.6e \n',iter,V_gate,norm( (dV)/length(dV) ));
    if iter > max_iter
        fprintf('Did not converge even after %d iterations\n',iter);
        check = 0;
    end
    
    if norm((dV)/length(dV)) < tolerance
        fprintf('Converged after %d iterations \n\n',iter);
        check = 0;
    end
    
end

%QQ(i)=1.6e-19*trapz(x,rho)/(Q*0.0529e-7*0.0529e-7);
Ec_new=Ec-Q*V0;
Ev_new=Ev-Q*V0;

if(OutGenerate.control.Solvers.SP)
    Eigen_val_C = OutSchr.Eigen_val_C;
    Eigen_val_V = OutSchr.Eigen_val_V;
    Psi2_C      = OutSchr.Psi2_C;
    Psi2_V      = OutSchr.Psi2_V;
end
if(OutGenerate.control.Solvers.SP)
    variables_name={'V_gate', 'V0','rho', 'n_charge', 'p_charge',...
        'Ec_new','Ev_new','OutGenerate','OutSchr','rho1'};
else
    variables_name={'V_gate', 'V0','rho', 'n_charge', 'p_charge',...
        'Ec_new','Ev_new','OutGenerate','rho1'};
end
if ~exist(fullfile(output_dir,'mat_files'),'dir')
    mkdir(fullfile(output_dir,'mat_files'))
end
filename = fullfile(output_dir,'mat_files',[file_prefix num2str(V_gate) '.mat']);
save(filename,variables_name{:});

% temporarily storing carrier densities for HFCV calculation
% use need modifications for IFCV (Intermidiate Frequency CV)
n_charge1=n_charge;
p_charge1=p_charge;


    


OutGenerate.control.Output_directory=out_dir{1};

end