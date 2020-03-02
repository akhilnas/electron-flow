%% GENERATE
% Generates parameters for simulation from input file
%% Update 
% - Delta_Ec is now obtained from electron affinity
% - Delta_Ec is now an unrecognized keyword

% - Update on 12th Aug 2014: 
%   OutGenerate.Ec, OutGenerate.Ev, OutGenerate.Ei nad OutGenerate.Ef are
%   now obtained by solving poisson only for equilibrium. Must when gate
%   stack has heterostructures. Equilibrium is solved using Neumann Neumann
%   boundary conditions. 

%% Function dfinition
function OutGenerate = generate (input_file)

%% Getting input file

if isequal(input_file,0)
    err = MException('MATLAB:InvalidFileFid','error: Unable to access the file');
    throw(err);
end
%% parser 

OutParser = parser(input_file,pwd);

%% OutParser: 

CheckOutParser(OutParser);

fprintf('Generating the variables...    \n')

var_surface   = OutParser.surface;
var_substrate = OutParser.substrate;
var_control   = OutParser.control;
var_layer     = OutParser.layer ;

%-------------------------------------------------------------------------
%-------------Band alignment generation -----------------------------------
%-------------needs better way of coding here -----------------------------

% Modifying Delta_Ec. Delta_Ec now obtained from electron affinity.
% removing Delta_Ec keyword. Delta_Ec is now unrecognized keyword.

for i=1:length(var_layer)-1
    var_layer(i).Delta_Ec=var_layer(i+1).layer_material.E_affinity-var_layer(i).layer_material.E_affinity;
end
var_layer(end).Delta_Ec=0;

Ec_temp = var_surface.barrier;
for i=1:size(var_layer,2)
    var_layer(i).Ec = Ec_temp;
    Ec_temp =  Ec_temp - var_layer(i).Delta_Ec;
end

%%   Generating the nodes and setting values to their properties

node_list = [];
node_list_temp = {};

for i=1:size(var_layer,2)
    
    temp_mat=[];
    for j=1:size(var_layer(i).sublayers,2)
        sl_start     = var_layer(i).Pos_begin + var_layer(i).sublayers(j).Pos_begin;
        sl_end       = var_layer(i).Pos_begin + var_layer(i).sublayers(j).Pos_end;
        temp_mat     = [temp_mat;[sl_start sl_end]];
        grid_spacing = var_layer(i).sublayers(j).grid_spacing;
        thickness    = var_layer(i).sublayers(j).thickness;
        No_nodes     = round(thickness/grid_spacing)+1;
        
        if (No_nodes > 0)
            temp_node(No_nodes) = node;
            
            for k=1:No_nodes
                temp_node(k).x           = sl_start + (k-1)*grid_spacing;
                temp_node(k).Nd          =  var_layer(i).sublayers(j).Nd;
                temp_node(k).Na          =  var_layer(i).sublayers(j).Na;
                temp_node(k).Nbt         =  var_layer(i).sublayers(j).Nbt;
%                 temp_node(k).temperature = var_layer(i).temperature;
                temp_node(k).layer       = i;
                temp_node(k).sublayer    = j;
                temp_node(k).epsilon     = var_layer(i).layer_material.epsilon;
                temp_node(k).mn_eff      = var_layer(i).layer_material.mn_eff;
                temp_node(k).mp_eff      = var_layer(i).layer_material.mp_eff;                
                temp_node(k).Eg          = var_layer(i).layer_material.Eg;
                temp_node(k).E_affinity  = var_layer(i).layer_material.E_affinity;
                temp_node(k).eta_C       = var_layer(i).layer_material.eta_C;
                temp_node(k).eta_V       = var_layer(i).layer_material.eta_V;
                temp_node(k).Ec          = var_layer(i).Ec;
                temp_node(k).Nc          = var_layer(i).layer_material.Nc;
                temp_node(k).Nv          = var_layer(i).layer_material.Nv;
            end
            
            node_list_temp{1+size(node_list_temp,1),1} = temp_node;
            clear temp_node;
        end
    end
    
    l_begin      = var_layer(i).Pos_begin;
    l_end        = var_layer(i).Pos_end;
    grid_spacing = var_layer(i).grid_spacing;
    
    for j=1:size(temp_mat,1)+1
        if (size(temp_mat,1)==0)
            thickness = abs(l_end-l_begin);
            No_nodes  = round(thickness/grid_spacing)+1;
            if(No_nodes > 0)
                temp_node(No_nodes) = node;
                
                for k=1:No_nodes
                    temp_node(k).x           = l_begin + (k-1)*grid_spacing;
                    temp_node(k).Nd          = var_layer(i).Nd;
                    temp_node(k).Na          = var_layer(i).Na;
                    temp_node(k).Nbt         = var_layer(i).Nbt;
%                     temp_node(k).temperature = var_layer(i).temperature;
                    temp_node(k).layer       = i;
                    temp_node(k).epsilon     = var_layer(i).layer_material.epsilon;
                    temp_node(k).mn_eff      = var_layer(i).layer_material.mn_eff;
                    temp_node(k).mp_eff      = var_layer(i).layer_material.mp_eff;                    
                    temp_node(k).Eg          = var_layer(i).layer_material.Eg;
                    temp_node(k).E_affinity  = var_layer(i).layer_material.E_affinity;
                    temp_node(k).eta_C       = var_layer(i).layer_material.eta_C;
                    temp_node(k).eta_V       = var_layer(i).layer_material.eta_V;
                    temp_node(k).Ec          = var_layer(i).Ec;
                    temp_node(k).Nc          = var_layer(i).layer_material.Nc;
                    temp_node(k).Nv          = var_layer(i).layer_material.Nv;
                end
                node_list_temp{1+size(node_list_temp,1),1} = temp_node;
                clear temp_node;
            end
            break;
        end
        
        if (sum(temp_mat(:,2)-temp_mat(:,1)) == var_layer(i).thickness)
            break;
        end
        
        if (j==1)
            
            thickness = abs(temp_mat(j,1)-l_begin);
            No_nodes  = round(thickness/grid_spacing);
        elseif (j==size(temp_mat,1)+1)
            
            thickness = abs(temp_mat(j-1,2)-l_end);
            No_nodes  = round(thickness/grid_spacing);
        else
            
            thickness = abs(temp_mat(j,1)-temp_mat(j-1,2));
            No_nodes  = round(thickness/grid_spacing);
        end
        
        if(No_nodes > 0)
            temp_node(No_nodes) = node;
            
            for k=1:No_nodes
                if (j==1)
                    temp_node(k).x = l_begin + (k-1)*grid_spacing;
                else
                    temp_node(k).x = temp_mat(j-1,2) + (k-1)*grid_spacing;
                end
                temp_node(k).Nd          = var_layer(i).Nd;
                temp_node(k).Na          = var_layer(i).Na;
                temp_node(k).Nbt         = var_layer(i).Nbt;
%                 temp_node(k).temperature = var_layer(i).temperature;
                temp_node(k).layer       = i;
                temp_node(k).epsilon     = var_layer(i).layer_material.epsilon;
                temp_node(k).mn_eff      = var_layer(i).layer_material.mn_eff;
                temp_node(k).mp_eff      = var_layer(i).layer_material.mp_eff;                
                temp_node(k).Eg          = var_layer(i).layer_material.Eg;
                temp_node(k).E_affinity  = var_layer(i).layer_material.E_affinity;
                temp_node(k).eta_C       = var_layer(i).layer_material.eta_C;
                temp_node(k).eta_V       = var_layer(i).layer_material.eta_V;
                temp_node(k).Ec          = var_layer(i).Ec;
                temp_node(k).Nc          = var_layer(i).layer_material.Nc;
                temp_node(k).Nv          = var_layer(i).layer_material.Nv;
            end
            node_list_temp{1+size(node_list_temp,1),1} = temp_node;
            clear temp_node;
        end
    end
    
end

sort_mat = [];
for i=1:size(node_list_temp,1)
    sort_mat = [sort_mat;[i node_list_temp{i}(1).x]];
end
sort_mat = sortrows(sort_mat,2);

for i=1:size(node_list_temp,1)
    node_list = [node_list node_list_temp{sort_mat(i,1)}] ;
end

temp=[];
for i=1:length(node_list)-1
    if ( (node_list(i).x)==(node_list(i+1).x) )
        node_list(i)=node_list(i+1);
        temp = [temp i+1];
        
    end
end
node_list(temp) = [];


%%   Generating variable vectors from nodes

global_variables;

% Defines all fundamental constants
Fundamental_Constants;

No_of_nodes = length(node_list);

OutGenerate.x              = zeros(1,No_of_nodes);
OutGenerate.epsilon_nodes  = zeros(1,No_of_nodes);
OutGenerate.mn_eff         = zeros(1,No_of_nodes);
OutGenerate.mp_eff         = zeros(1,No_of_nodes);
OutGenerate.Ec             = zeros(1,No_of_nodes);
OutGenerate.Ev             = zeros(1,No_of_nodes);
OutGenerate.Na             = zeros(1,No_of_nodes);
OutGenerate.Nd             = zeros(1,No_of_nodes);
OutGenerate.Nbt            = zeros(1,No_of_nodes);
OutGenerate.Nc             = zeros(1,No_of_nodes);
OutGenerate.Nv             = zeros(1,No_of_nodes);
OutGenerate.eta_C          = zeros(1,No_of_nodes);
OutGenerate.eta_V          = zeros(1,No_of_nodes);
OutGenerate.Eg             = zeros(1,No_of_nodes);
OutGenerate.Ec             = zeros(1,No_of_nodes);
OutGenerate.Ev             = zeros(1,No_of_nodes);

for i=1:No_of_nodes
    OutGenerate.x(i)             =18.9* node_list(i).x;
    OutGenerate.epsilon_nodes(i) = node_list(i).epsilon;
    OutGenerate.mn_eff(i)        = node_list(i).mn_eff *m_e;
    OutGenerate.mp_eff(i)        = node_list(i).mp_eff *m_e;
    OutGenerate.Ec(i)            = node_list(i).Ec /13.6;
    OutGenerate.Ev(i)            = node_list(i).Ev /13.6;
    OutGenerate.Na(i)            = node_list(i).Na *Bohr_radius^3;
    OutGenerate.Nd(i)            = node_list(i).Nd *Bohr_radius^3;
    OutGenerate.Nbt(i)           = node_list(i).Nbt *Bohr_radius^3;
    OutGenerate.eta_C(i)         = node_list(i).eta_C;
    OutGenerate.eta_V(i)         = node_list(i).eta_V;
    OutGenerate.Eg(i)            = node_list(i).Eg /13.6;    
    OutGenerate.Nc(i)            = node_list(i).Nc *Bohr_radius^3*((OutParser.Temperature/300)^(3/2));
    OutGenerate.Nv(i)            = node_list(i).Nv *Bohr_radius^3*((OutParser.Temperature/300)^(3/2));
%     OutGenerate.Nc1(i)            = (node_list(i).eta_C)*[2* (2*pi*(node_list(i).mn_eff *m_e) *Kb*(( KbT/ 300) * node_list(i).temperature)/h^2 ).^(3/2)];
%     OutGenerate.Nv1(i)            = 2* (2*pi*(node_list(i).mp_eff *m_e) *Kb*(( KbT/ 300) * node_list(i).temperature)/h^2 ).^(3/2);
    
end

% OutGenerate.Tempr = ( KbT/ 300) * node_list(1).temperature;
OutGenerate.Tempr = ( KbT/ 300) * OutParser.Temperature;
for i=1:No_of_nodes
    OutGenerate.Nc1(i)            = (node_list(i).eta_C)*2* (2*pi*(node_list(i).mn_eff *m_e) *Kb*( OutGenerate.Tempr)/h^2 ).^(3/2);
    OutGenerate.Nv1(i)            = 2* (2*pi*(node_list(i).mp_eff *m_e) *Kb*( OutGenerate.Tempr)/h^2 ).^(3/2);
end

OutGenerate.layer     = var_layer;
OutGenerate.control   = var_control;
OutGenerate.surface   = var_surface;
OutGenerate.substrate = var_substrate;

%%  Generating dependent parameters
mn_eff   = OutGenerate.mn_eff;
mp_eff   = OutGenerate.mp_eff;
Tempr = OutGenerate.Tempr;
eta_C=OutGenerate.eta_C;
eta_V=OutGenerate.eta_V;

%for i=1:No_of_nodes
%Nc = (eta_C(i))*[2* (2*pi*mn_eff *Kb*Tempr/h^2 ).^(3/2)] ; % Nc per Bohr_rad^3, accounting for degeneracy
%end
%Nv= 2* (2*pi*mp_eff *Kb*Tempr/h^2 ).^(3/2);  % Nc per Bohr_rad^3


Nc    = OutGenerate.Nc;
Nv    = OutGenerate.Nv;
Nd    = OutGenerate.Nd;
Na    = OutGenerate.Na;
Eg    = OutGenerate.Eg;
Ec    = OutGenerate.Ec;


ni = sqrt(Nc.*Nv).*exp(-Eg/(2*Kb*Tempr));  % Intrinsic carrier density per Bohr_rad^3

Ei = OutGenerate.Ec + Kb*Tempr*log(ni./Nc);

%Calculating Ef
Ef = zeros(1,No_of_nodes);
a = Nd > Na;
Ef(a) = Ei(a) + Kb*Tempr*log((Nd(a)-Na(a))./ni(a));
a = Nd < Na;
Ef(a) = Ei(a) - Kb*Tempr*log((Na(a)-Nd(a))./ni(a));
a = Nd == Na;
Ef(a) = Ei(a);
% Ef =  Ef(end) * ones(1,No_of_nodes);
Ef_=Ef(~isnan(Ef));
OutGenerate.Ec = OutGenerate.Ec - Ef_(end);
OutGenerate.Ev = OutGenerate.Ev - Ef_(end);

Ei = Ei - Ef(end);
Ef=Ef-Ef(end);

epsilon = zeros(1,No_of_nodes-1);
for i=1:No_of_nodes-1
    epsilon(i) = 0.5*(OutGenerate.epsilon_nodes(i)+OutGenerate.epsilon_nodes(i+1));
end

OutGenerate.ni = ni;
OutGenerate.Ef = Ef;
OutGenerate.Ei = Ei;
OutGenerate.epsilon = epsilon;

% Saving "just joining" band diagram with '_old' suffix
OutGenerate.Ec_old=OutGenerate.Ec;
OutGenerate.Ev_old=OutGenerate.Ev;
OutGenerate.Ei_old=OutGenerate.Ei;
OutGenerate.Ef_old=OutGenerate.Ef;

OutGenerate.Ef = zeros(1,length(Ef));
H=getH(OutGenerate);
OutGenerate=getEquilibrium(H,OutGenerate);

%%  Calculating Phi_ms
Phi_ms = 0;
% if var_substrate.Zero_Slope && var_surface.Schottky
%     Phi_ms = var_surface.barrier - 19.2*( OutGenerate.Ec(end) - OutGenerate.Ef(end));
% end

OutGenerate.Phi_ms = Phi_ms;

%%  Setting Schrodinger properties

if var_control.Solvers.SP
    % finding the closest node point to SchodingerStart and SchrodingerEnd
    [III , SchrStart]  = min(abs((0.0529*OutGenerate.x - var_control.Schodinger_Poisson.SchrodingerStart)));
    [III , SchrStop]   = min(abs((0.0529*OutGenerate.x - var_control.Schodinger_Poisson.SchrodingerStop )));
    if (SchrStop - SchrStart) < 5
        error('Please choose the range for Schrodinger solver in such a way that there are at least 5 grid points')
    end
    
    OutGenerate.SchrStart = SchrStart;
    OutGenerate.SchrStop  = SchrStop;    
end

fprintf('Done \n\n')

end



