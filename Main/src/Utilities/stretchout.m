clc;
clear all;
close all;

%% Importing measured CV data;
[fn pn]=uigetfile('*.csv','please select C-V file');
file=fullfile(pn,fn);
CVm=csvread(file,1,0);
Vm=CVm(:,1);
Cm=CVm(:,2);

% making all values in Cm distinct
    for i=2:length(Cm)
        %     if Cm(i)==Cm(i-1);
        %         Cm(i)
        %         Cm(i-1)
        %         i
        %          Cm(i)=Cm(i)+0.00001;
        %     end
        
        zzz=find(not((Cm)-Cm(i)));
        
        if length(zzz)>1
            for j=1:length(zzz)
                Cm(zzz(j))=Cm(zzz(j))+zzz(j)*1e-7;
            end
        end
        
    end
Cm=Cm*1e-6;

%% Importing Ideal CV
Ideal_Dir=uigetdir(pwd);
Tox=input('please provide dielectric thickness in nm\n','s');
Tox=eval(Tox);
Tox=Tox*18.9;
Cox=3.9*8.85e-14/(Tox*1e-7/18.9); % use only when ideal simulations are with SiO2 only dielectric
[OutGenerate HFCV]=Ideal_HFCV(Ideal_Dir,Tox);
Vg=HFCV.Vg;
C=HFCV.C;
psi_s=HFCV.psi_s;
Vpsi=HFCV.psi_s/13.6;

[forget sur_loc_index]=min(abs(OutGenerate.x-Tox));

Ei_s=OutGenerate.Ei(sur_loc_index)-Vpsi;

%% equaling data points in ideal and measured
Cm1=zeros(length(C),1);
Vm1=zeros(length(Vg),1);

for i = 1:length(C)
    Vm1(i)=interp1(Cm,Vm,C(i));
end

Cm1=C;

%% Calculating Dit

Cox=HFCV.Cox;

smooth(Vm1);

D1=diff(Vg)./diff(psi_s);
D2=diff(Vm1)./diff(psi_s);
Cit=Cox*(D2-D1);
Dit=Cit/1.6e-19;

Ef_Ev=(OutGenerate.Ef(sur_loc_index)-OutGenerate.Ev(sur_loc_index))*13.6;
E=psi_s(1:end-1)+Ef_Ev;


%% Plotting

ZBB=interp1(psi_s,Vg,0);
if (isnan(ZBB))
    ZBB=-1.5;
end
CFB=interp1(Vg,C,ZBB);
VFB=interp1(Cm,Vm,CFB);
Vm=Vm-VFB+ZBB;

subplot(1,2,1);
plot(Vg,C,'b','linewidth',2);
hold on;
%     plot(Vm,Cm*1e-6,'k');
plot(Vm,Cm,'k','linewidth',2);
hold off;
set(gca,'fontsize',20)
xlabel('gate bias (V)');
ylabel('capacitance (F cm^{-2})');
xlim([-2.5 2]);
grid on;
set(gca,'gridlinestyle','--');
legend('Ideal HFCV','1 MHz');

subplot(1,2,2);
plot(E,Dit,'sb');
set(gca,'fontsize',20)
xlabel('E-Ev (eV)');
ylabel('D_{it} (eV^{-1} cm^{-2})');
xlim([-0.4 0.6]);
grid on;
set(gca,'gridlinestyle','--');

%%


% Vmg=interp1(Ei_s,Vg,0); % midgap voltage in ideal CV
% Cmg=interp1(Vg,C,Vmg); % midgap capacitance in ideal and measured CV
% 
% Vmg_m=interp1(Cm,Vm,Cmg); % midgap voltage in measured CV
% 
% Vfix=Vmg-Vmg_m;
% 
% 
% plot(Vg,C,Vm+Vfix,Cm)