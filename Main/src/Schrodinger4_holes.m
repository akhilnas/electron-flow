function  OutSchr = Schrodinger4_holes(OutGenerate,Ec_new,Ev_new)

UseSchr.V = 1;
N_eig=15;
global hbar 

persistent  Hv  Hv_diag Psi_V Psi2_V L invL IsFirst  Ldiag d

x = OutGenerate.x;
mp_eff = OutGenerate.mp_eff;
SchrStart = OutGenerate.SchrStart;
SchrStop = OutGenerate.SchrStop;
Alpha = 2;
if isempty(IsFirst)
    IsFirst = 0;
end
%%  Initialising Variables

if IsFirst == 0
    
    No_of_nodes_schr = SchrStop - SchrStart + 1;
    size_of_matrix = No_of_nodes_schr-2;
    
    Hupper = zeros(1,size_of_matrix-1);
    Hlower = zeros(1,size_of_matrix-1);
    Ldiag  = zeros(1,size_of_matrix);
    Hv_diag = zeros(1,size_of_matrix);
    
    mp_eff_region = mp_eff(1:end-1);
    
    Psi_V  = zeros(size_of_matrix+2,size_of_matrix);
    Psi2_V = zeros(size_of_matrix+2,size_of_matrix);
    
    d = diff(x);
end

%%  Schrodinger Equation for Valence band

if (IsFirst == 0)
    
    for n=2:No_of_nodes_schr-1
        Ldiag(n-1)=sqrt( (d(n+SchrStart-1)+d(n+SchrStart-2))/2  ) ;
    end
    L = diag(Ldiag,0);
    invL = diag(1./Ldiag,0); % faster and accurate way of finding inv of a diagonal matrix rather than using INV();
    
    for n=2:No_of_nodes_schr-1
        temp1=1/(mp_eff_region(n+SchrStart-2)*d(n+SchrStart-2));
        temp2=1/(mp_eff_region(n+SchrStart-1)*d(n+SchrStart-1));
        
        Hv_diag(n-1)=Alpha*(hbar^2/(d(n+SchrStart-1)+d(n+SchrStart-2))) *(temp1+temp2);
        
        if n<=No_of_nodes_schr-2
            Hupper(n-1)=-Alpha*(hbar^2/(d(n+SchrStart-1)+d(n+SchrStart-2)))*temp2;
        end
        
        if n>=3
            Hlower(n-2)=-Alpha*(hbar^2/(d(n+SchrStart-1)+d(n+SchrStart-2)))*temp1;
        end
    end
    
    Hv = diag(Hv_diag,0)+diag(Hupper,1)+diag(Hlower,-1);
    Hv = L*Hv*invL;
    Hv(2:size(Hv,1)+1:end) = diag(Hv,1);
    % Equating the lower diagonal to the upper
    % diagonal to expunge the small neumerical error
    
end

if UseSchr.V
    
    Hv(1:size(Hv,1)+1:end) = Hv_diag - Ev_new(1+SchrStart:SchrStop-1);
    
     if any(isnan(Ev_new))
        msgStr = 'error:Convergence failed';
        err = MException('MATLAB:ConvergenceFailed',msgStr);
        throw(err);
     end
    

    
    HHV=sparse(Hv);
    options.disp=0;
    [Psi_V(2:end-1,1:N_eig) Eigen_val_matrix] = eigs(HHV,N_eig,min(-Ev_new),options);
    
    for i=1:size(Psi_V(2:end-1,:),2)
        Psi2_V(2:end-1,i) = ( conj(Psi_V(2:end-1,i)./Ldiag') ).* (Psi_V(2:end-1,i)./Ldiag');
    end
    
    Eigen_val_V = -diag(Eigen_val_matrix);  % in Rydbergs


end

IsFirst = IsFirst + 1;

OutSchr.Psi2_V = Psi2_V;
OutSchr.Eigen_val_V = Eigen_val_V;


%% Plotting Wavefunction
% Conduction Band

DrawPlot = 0;
PlotWavefn = 1;
if (DrawPlot && PlotWavefn)
    
    subplot(3,2,1)
    for i=1:3
        temp=Ec_new*13.6;
        %      temp=temp/range(temp);
        yy1=plot(x*0.0529,temp);
        hold on
        temp1=Psi2_C(:,i);
        temp1=temp1/(max(temp1));
        yy2=plot(x(SchrStart:SchrStop)'*0.0529,temp1+temp(1));
        hold on
    end
    refreshdata(yy1,'caller')
    refreshdata(yy2,'caller')
    drawnow,pause(.001)
    
    plot(x*0.0529,Ec_new*13.6);
    hold off
    
    %Valence Band
    subplot(3,2,2)
    for i=1:3
        
        temp=Ev_new*13.6;
        %      temp=temp/range(temp);
        yy3=plot(x*0.0529,temp);
        hold on
        temp1=Psi2_V(:,i);
        temp1=-temp1/(max(temp1));
        yy4=plot(x(SchrStart:SchrStop)'*0.0529,temp1+temp(1));
        hold on
    end
    refreshdata(yy3,'caller')
    refreshdata(yy4,'caller')
    drawnow,pause(.01)
    
    plot(x*0.0529,Ec_new*13.6 );
    hold off
    
end


end