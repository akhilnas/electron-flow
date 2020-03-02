function OutGenerate=getEquilibrium(H,OutGenerate)

fprintf('Obtaining bands in equilibrium \n(also considering Phi_ms=0 for this analysis)\n');
global Q eps0 Kb;
Vt=Kb*OutGenerate.Tempr;

H=sparse(H);

max_iter         = 1000;
tolerance        = OutGenerate.control.Schodinger_Poisson.tolerance;
output_dir       = OutGenerate.control.Output_directory;
file_prefix      = OutGenerate.control.file_prefix;
NoMix            = 0;
alpha            = 1;
Ec               = OutGenerate.Ec;
Ev               = OutGenerate.Ev;
No_of_nodes      = length(OutGenerate.x);
x                = OutGenerate.x;
epsilon          = OutGenerate.epsilon;
Ef               = OutGenerate.Ef;
Nc               = OutGenerate.Nc;
Nv               = OutGenerate.Nv;
Tempr            = OutGenerate.Tempr;

% Vt               = OutGenerate.Tempr; % Boltzmann constant is Kb=1 in Rydberg units
N=length(x);


%Neumann Neumann Boundary conditions
tmp=H(1,1);
H(1,:)=0;
H(1,1)=tmp/2;
H(1,2)=-tmp/2;
H(2,1)=-tmp/2;

tmp=H(end,end);
H(end,:)=0;
H(end,end)=tmp/2;
H(end,end-1)=-tmp/2;
H(end-1,end)=-tmp/2;

V0 = zeros(1,No_of_nodes);

iter=0;
check=1;
V0(1)=0;
V0(end)=0;

while(check)
    iter=iter+1;
    
    n_charge = Nc.*exp((Ef-(Ec-Q*V0))/(Kb*Tempr) );
    p_charge = Nv.*exp(((Ev-Q*V0)-Ef)/(Kb*Tempr) );
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
    
    if NoMix==0
        V0=V0+alpha*dV';
    else
        V1=V0(1);
        V01 = Anderson_mixing(V0,V0+dV',iter,NoMix,alpha);
        V0=V01;
        V0(1)=V1;
    end
    
    fprintf('Iteration = %d \t Bias = equilibrium \t Error in V = %.6e \n',iter,norm( (dV)/length(dV) ));
    
    if iter > max_iter
        fprintf('Did not converge even after %d iterations\n',iter);
        check = 0;
    end
    
    if norm((dV)/length(dV)) < tolerance
        fprintf('Converged after %d iterations \n\n',iter);
        check = 0;
    end
    
end

QQ=1.6e-19*trapz(x,rho)/(Q*0.0529e-7*0.0529e-7);
Ec_new=Ec-Q*V0;
Ev_new=Ev-Q*V0;
V_gate=0;
OutGenerate.Ec=Ec_new;
OutGenerate.Ev=Ev_new;
OutGenerate.Ei=OutGenerate.Ei-Q*V0;

variables_name={'V_gate', 'V0','rho', 'n_charge', 'p_charge',...
    'Ec_new','Ev_new','OutGenerate'};

if ~exist(fullfile(output_dir,'mat_files','Equilibrium'),'dir')
    mkdir(fullfile(output_dir,'mat_files','Equilibrium'))
end


filename = fullfile(output_dir,'mat_files','Equilibrium',[file_prefix '_equilibrium.mat']);
fprintf('Saving results of equlibrium analysis in %s \n',filename);
save(filename,variables_name{:});


end