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

function scs_new(H,OutGenerate,V_ramp)


global Q eps0 Kb;
Vt=Kb*OutGenerate.Tempr;

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
        % Calculation for holes
        OutSchr_common = Schrodinger4_holes(OutGenerate,(Ec-Q*V0),(Ev-Q*V0));
        % Calculation for electrons
        mn_t      = OutGenerate.layer(1,2).layer_material.mn_t*ones(1,No_of_nodes);
        OutSchr_t = Schrodinger4_electrons(OutGenerate,(Ec-Q*V0),(Ev-Q*V0),mn_t);
        OutSchr_t.Psi2_V = OutSchr_common.Psi2_V;
        OutSchr_t.Eigen_val_V = OutSchr_common.Eigen_val_V;
        OutSchr_t.direction = 't';
        [rho_t, n_charge_t, p_charge, rho1] = charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0),OutSchr_t,0);
        mn_l      = OutGenerate.layer.layer_material.mn_l*ones(1,No_of_nodes);
        OutSchr_l = Schrodinger4_electrons(OutGenerate,(Ec-Q*V0),(Ev-Q*V0),mn_l);
        OutSchr_l.Psi2_V = OutSchr_common.Psi2_V;
        OutSchr_l.Eigen_val_V = OutSchr_common.Eigen_val_V;
        OutSchr_l.direction = 'l';
        [rho_l, n_charge_l, p_charge, rho1] = charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0),OutSchr_l,0);
        n_charge = n_charge_t + n_charge_l;
        rho = (rho_t + rho_l - n_charge)/2 + n_charge;
    else
        [rho n_charge p_charge rho1]=charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0),0,1);
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
    Eigen_val_C_t = OutSchr_t.Eigen_val_C;
    Eigen_val_C_l = OutSchr_l.Eigen_val_C;
    Eigen_val_V = OutSchr_common.Eigen_val_V;
    Psi2_C_t      = OutSchr_t.Psi2_C;
    Psi2_C_l      = OutSchr_l.Psi2_C;
    Psi2_V      = OutSchr_common.Psi2_V;
end
if(OutGenerate.control.Solvers.SP)
    variables_name={'V_gate', 'V0','rho', 'n_charge', 'p_charge',...
        'Ec_new','Ev_new','OutGenerate','OutSchr_t','OutSchr_l','rho1'};
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

% Random Plots
% figure(1);
% plot(x*0.0529,n_charge/((5.29*10^-9)^3));
% xlabel('Distance(nm)');
% ylabel('Electron Density(cm-3)');
% legend('Single Effective Mass','Two Mass Solver');
% hold all;
% figure(2);
% plot(x*0.0529,V0*19.2);
% xlabel('Distance(nm)');
% ylabel('Potential(V)');
% legend('Single Effective Mass','Two Mass Solver');
% hold all;
    


OutGenerate.control.Output_directory=out_dir{1};

end