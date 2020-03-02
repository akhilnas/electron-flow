global_variables;

Tempr=KbT/300 *TemprK(1);  % Temperature in atomic units

No_of_regions = No_of_nodes-1;
d = zeros(1,No_of_regions);

for i=1:No_of_nodes-1
    d(i) = x(i+1)-x(i);
end

% mn_eff_node=zeros(1,No_of_nodes);
% mp_eff_node=zeros(1,No_of_nodes);
% for i=2:No_of_nodes-1
%     mn_eff_node(i)=0.5*(mn_eff(i-1)+mn_eff(i));
%     mp_eff_node(i)=0.5*(mp_eff(i-1)+mp_eff(i));
% end
% mn_eff_node(1)=mn_eff(1);mn_eff_node(No_of_nodes)=mn_eff(No_of_nodes-1);
% mp_eff_node(1)=mp_eff(1);mp_eff_node(No_of_nodes)=mp_eff(No_of_nodes-1);

epsilon = zeros(1,No_of_nodes-1);
for i = 1:No_of_nodes-1
    epsilon(i)=0.5*(epsilon_nodes(i)+epsilon_nodes(i+1));
end

Nc = 2* (2*pi*mn_eff *Kb*Tempr/h^2 ).^(3/2);  % Nc per Bohr_rad^3
Nv = 2* (2*pi*mp_eff *Kb*Tempr/h^2 ).^(3/2);  % Nc per Bohr_rad^3

ni = sqrt(Nc.*Nv).*exp(-Eg/(2*KbT));  % Intrinsic carrier density per Bohr_rad^3

Ei=Ec_nodes+Kb*Tempr*log(ni./Nc);


for i=1:No_of_nodes
    if Nd(i)-Na(i) > 0
        Ef=Ei+KbT*log((Nd(i)-Na(i))./ni(i));
    elseif Na(i)-Nd(i) >0
        Ef=Ei-KbT*log((Na(i)-Nd(i))./ni(i));
    elseif Na(i)-Nd(i)==0
        Ef(i)=Ei(i);
    end
end

% This is a temporary code


Ef =  Ef(end) * ones(1,No_of_nodes);

Ec_nodes = Ec_nodes -Ef;
Ev_nodes = Ev_nodes -Ef;
Ei = Ei-Ef;
Ef = zeros(1,length(Ef));

% figure;
% plot((x-x(Oxide_nodes(1)-1))*0.0529,Ec_nodes*13.6);hold on
% plot((x-x(Oxide_nodes(1)-1))*0.0529,Ev_nodes*13.6);hold on
% x_metal = x(end) + 18.9*[0:0.2:5];
% plot((x_metal-x(Oxide_nodes(1)-1))*0.0529,0.17*ones(size(x_metal)));hold off
% 
% hold on;
% plot((x-x(Oxide_nodes(1)-1))*0.0529,Ef*13.6)
% 1;

 
%%%%%%%%%%%%%%%% JOYS-DIXON APPROXIMATION FOR FERMI LEVEL %%%%%%%%%%%%%%%%%

%Ef=Ev_nodes-KbT*( log(Na./Nv) + (1/sqrt(8))* Na./Nv );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( exist('SchrodingerStart','var') && exist('SchrodingerEnd','var') ...
       && (SchrodingerStart < SchrodingerEnd) ) 
    UseSchrodinger = 1;
% finding the closest node point to SchodingerStart and SchrodingerEnd 
    [notReqrd, SchrStart] = min(abs((0.0529*x - SchrodingerStart)));
    [notReqrd, SchrEnd]   = min(abs((0.0529*x - SchrodingerEnd  )));
    if (SchrEnd-SchrStart) < 5
        error('Please choose the range for Schrodinger solver in such a way that there are at least 5 grid points')
    end
end 

% temporary
UseSchrodinger = 1;
SchrStart = 1;
SchrEnd   = No_of_nodes;




