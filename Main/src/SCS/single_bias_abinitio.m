function single_bias_abinitio()
% Density-Gradient Solver with Fitting Based on Ab-initio Data
% Started on 01 Feb 2018 by Akhil
% Maintained by : Akhil (eakhil711@gmail.com)

% All Input Files of Different Widths
files_SP = {'stack_files/SP/3_2nm.txt' , 'stack_files/SP/3_9nm.txt' , 'stack_files/SP/5_4nm.txt' , 'stack_files/SP/8_7nm.txt' ,...
    'stack_files/SP/10_3nm.txt' };
files_DG = {'stack_files/DG/3_2nm.txt' , 'stack_files/DG/3_9nm.txt' , 'stack_files/DG/5_4nm.txt' , 'stack_files/DG/8_7nm.txt' ,...
    'stack_files/DG/10_3nm.txt' };

% Selection of n-type(n) or p-type(p) simulation
Simulation_type = input('Enter [n] for n-type simulation or [p] for p-type simulation\n','s');

% Selection of Simulation Procedure
% Optimization of Fitting Factors to Ab-initio data(optimization)
% Sweep of Fitting factor over various values(sweep)
Simulation_procedure = input('Enter [optimization] for Optimization of Fitting factors to Ab-initio data or enter [sweep] for a sweep of fitting factor values\n','s');

% Choice on running Schrodinger-Poisson solver
choice = input('Do you wish to run Schrodinger-Poisson. Reply with [Y] to run \n','s');


for fileloop = 1:length(files_SP)
    
    % Processing of Input File for Schrodinger-Poisson solver
    input_file_SP = files_SP{fileloop};
    OutGenerate=generate(input_file_SP); 
    H=getH(OutGenerate);
    H=sparse(H);    

    % Selection of Appropriate Ab-initio data for current width Simulation
    % 1 for 3.2 nm
    % 2 for 3.9 nm
    % 3 for 5.4 nm
    % 4 for 8.7 nm    
    % 5 for 10.3 nm
    wchoice = fileloop;    
    
    if Simulation_type == 'n'
        if wchoice == 4
            % Ab-initio data for 8.7 nm fin width (n-type)
            Extra_electron = [-0, -0.00020408, -0.00040816, -0.00061224, -0.00081633, -0.00102041, -0.00122449, -0.00142857,... 
                -0.00163265, -0.00183673, -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714, -0.00306122,... 
                -0.00326531, -0.00346939, -0.00367347, -0.00387755, -0.00408163, -0.00428571, -0.0044898 , -0.00469388,... 
                -0.00489796, -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837, -0.00612245, -0.00632653,... 
                -0.00653061, -0.00673469, -0.00693878, -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,... 
                -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959, -0.00918367, -0.00938776, -0.00959184,... 
                -0.00979592, -0.01];
            Gate_Bias = [0,  0.42212308,  0.44607946,  0.46255743,  0.47595953, 0.48765822,  0.49826737,  0.50811615,  0.5173989,...  
                0.52624297, 0.53473513,  0.5429375 ,  0.55089615,  0.55864624,  0.56621529, 0.57362521,  0.58089533,  0.58803734,... 
                0.595065  ,  0.6019887 , 0.6088173 ,  0.61555846,  0.62221998,  0.62880525,  0.63532063, 0.64177207,  0.64816184,...  
                0.65449354,  0.66077104,  0.66699734, 0.67317528,  0.67930715,  0.68539533,  0.69144191,  0.69744886, 0.70341794,...  
                0.70935083,  0.71524906,  0.72111408,  0.72694725, 0.73274979,  0.73852292,  0.74426772,  0.74998524,  0.75567646,...
                0.76134231,  0.76698366,  0.77260134,  0.77819613,  0.78376877];
        elseif wchoice == 3
            % Ab-initio Data for 5.46 nm fin width (n-type)
            Extra_electron = [-0, -0.00020408, -0.00040816, -0.00061224, -0.00081633, -0.00102041, -0.00122449, -0.00142857,... 
                -0.00163265, -0.00183673, -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714, -0.00306122,... 
                -0.00326531, -0.00346939, -0.00367347, -0.00387755, -0.00408163, -0.00428571, -0.0044898 , -0.00469388,... 
                -0.00489796, -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837, -0.00612245, -0.00632653,... 
                -0.00653061, -0.00673469, -0.00693878, -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,...
                -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959, -0.00918367, -0.00938776, -0.00959184,... 
                -0.00979592, -0.01];
            Gate_Bias = [0,  0.44933436,  0.47257258,  0.4883634 ,  0.5010991 , 0.5121562 , 0.52214786,  0.53140517,  0.54011577,...  
                0.54840932, 0.5563719 ,  0.5640628 ,  0.57153152,  0.57881053,  0.5859266 , 0.59290101,  0.59975093,  0.6064904 ,...  
                0.6131311 ,  0.6196828 , 0.62615379,  0.63255117,  0.63888159,  0.64514928,  0.65135941, 0.65751605,  0.66362283,...  
                0.66968298,  0.67569935,  0.68167453, 0.68761085,  0.69351043,  0.69937516,  0.7052068 ,  0.71100695, 0.71677663,...  
                0.72251815,  0.72823211,  0.73391979,  0.73958226, 0.7452205 ,  0.75083545,  0.75642763,  0.76199875,  0.76754836,...
                0.77307821,  0.77858859,  0.78408011,  0.78955337,  0.79500894];
        elseif wchoice == 1    
            % Ab-initio Data for 3.2 nm fin width (n-type)
            Extra_electron = [-0, -0.00020408, -0.00040816, -0.00061224, -0.00081633, -0.00102041, -0.00122449, -0.00142857,... 
                -0.00163265, -0.00183673, -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714, -0.00306122,... 
                -0.00326531, -0.00346939, -0.00367347, -0.00387755, -0.00408163, -0.00428571, -0.0044898 , -0.00469388,... 
                -0.00489796, -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837, -0.00612245, -0.00632653,... 
                -0.00653061, -0.00673469, -0.00693878, -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,...
                -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959, -0.00918367, -0.00938776, -0.00959184,... 
                -0.00979592, -0.01];
            Gate_Bias = [0,  0.48374189,  0.50645215,  0.52172187,  0.53394579, 0.54450006,  0.55399781,  0.5627661 ,  0.57100015,...  
                0.57882589, 0.5863293 ,  0.5935745 ,  0.60060139,  0.60744706,  0.61413807,  0.62069256,  0.6271357 ,  0.63347464,...  
                0.6397224 ,  0.64588923, 0.65198315,  0.65801115,  0.66397925,  0.66989268,  0.67575598,  0.68157314,  0.68734768,...  
                0.69308273,  0.69878106,  0.70444517, 0.71007729,  0.71567942,  0.7212534 ,  0.72680087,  0.73232334,  0.73782219,...  
                0.74329866,  0.74875391,  0.754189  ,  0.75960491, 0.76500254,  0.77038271,  0.77574622,  0.78109377,  0.78642604,...
                0.79174363,  0.79704715,  0.80233713,  0.80761407,  0.81287846];
        elseif wchoice == 2   
            % Ab-initio Data for 3.9 nm fin width (n-type)
            Extra_electron = [-0, -0.00020408, -0.00040816, -0.00061224, -0.00081633, -0.00102041, -0.00122449, -0.00142857,... 
                -0.00163265, -0.00183673, -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714, -0.00306122,... 
                -0.00326531, -0.00346939, -0.00367347, -0.00387755, -0.00408163, -0.00428571, -0.0044898 , -0.00469388,... 
                -0.00489796, -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837, -0.00612245, -0.00632653,... 
                -0.00653061, -0.00673469, -0.00693878, -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,...
                -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959, -0.00918367, -0.00938776, -0.00959184,... 
                -0.00979592, -0.01];
            Gate_Bias = [0,  0.47377065,  0.4966183 ,  0.51202585,  0.52438554, 0.53507386,  0.54470391,  0.55360277,  0.56196564,...  
                0.56991844, 0.57754713,  0.58491305,  0.59206145,  0.59902673,  0.6058356 ,  0.61250924,  0.61906674,  0.62551796,...  
                0.63187657,  0.63815258, 0.64435356,  0.65048682,  0.65655812,  0.66257344,  0.66853656,  0.67445171,  0.68032241,...  
                0.68615178,  0.69194239,  0.69769678, 0.7034177 ,  0.7091068 ,  0.71476594,  0.72039673,  0.72600069,  0.73157919,...  
                0.73713402,  0.74266526,  0.74817452,  0.75366228, 0.75913048,  0.76457938,  0.77000978,  0.77542242,  0.78081798,...
                0.78619707,  0.79156029,  0.79690819,  0.80224127,  0.80756004];
        elseif wchoice == 5    
            % Ab-initio Data for 10.3 nm fin width (n-type)
            Extra_electron = [-0, -0.00020408, -0.00040816, -0.00061224, -0.00081633, -0.00102041, -0.00122449, -0.00142857,... 
                -0.00163265, -0.00183673, -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714, -0.00306122,... 
                -0.00326531, -0.00346939, -0.00367347, -0.00387755, -0.00408163, -0.00428571, -0.0044898 , -0.00469388,... 
                -0.00489796, -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837, -0.00612245, -0.00632653,... 
                -0.00653061, -0.00673469, -0.00693878, -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,...
                -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959, -0.00918367, -0.00938776, -0.00959184,... 
                -0.00979592, -0.01 ];
            Gate_Bias = [0,  0.41619826,  0.44048321,  0.45727504,  0.47096818, 0.48294097,  0.49380834,  0.50390104,  0.51341307,...  
                0.52247307, 0.5311686 ,  0.53956254,  0.54770167,  0.55562186,  0.56335936,  0.57091603,  0.57832723,  0.58560369,...  
                0.59275826,  0.59980181, 0.60674939,  0.61359201,  0.62035398,  0.62703555,  0.63364227,  0.64017891,  0.64665301,...  
                0.6530625 ,  0.65941392,  0.6657105 , 0.67195519,  0.67815068,  0.68429947,  0.69040383,  0.69646588,  0.70248757,...  
                0.7084707 ,  0.71441698,  0.72032796,  0.72620511, 0.73204981,  0.73786334,  0.74364691,  0.74940165,  0.75512863,...
                0.76082885,  0.76650225,  0.77215106,  0.77777883,  0.78337977];
        end
    elseif Simulation_type == 'p'
        if wchoice == 1
            % Ab-initio Data for 3.2 nm fin width (p-type)
            Extra_hole = [0,  0.00020408,  0.00040816,  0.00061224,  0.00081633,  0.00102041,  0.00122449,  0.00142857,...  
                0.00163265,  0.00183673,  0.00204082,  0.0022449 ,  0.00244898,  0.00265306,  0.00285714,  0.00306122,...  
                0.00326531,  0.00346939,  0.00367347,  0.00387755,  0.00408163,  0.00428571,  0.0044898 ,  0.00469388,...  
                0.00489796,  0.00510204,  0.00530612,  0.0055102 ,  0.00571429,  0.00591837,  0.00612245,  0.00632653,...  
                0.00653061,  0.00673469,  0.00693878,  0.00714286,  0.00734694,  0.00755102,  0.0077551 ,  0.00795918,...
                0.00816327,  0.00836735,  0.00857143,  0.00877551,  0.00897959,  0.00918367,  0.00938776,  0.00959184,...  
                0.00979592,  0.01];
            Gate_Bias = [0, -0.45879401, -0.48130773, -0.49638348, -0.5084134 , -0.51877399, -0.52807836, -0.53665356,... 
                -0.54469483, -0.55232808, -0.5596393 , -0.56668981, -0.57352488, -0.58017891, -0.58667863, -0.59304522,... 
                -0.59929577, -0.60544424, -0.61150221, -0.61747938, -0.62338394, -0.62922291, -0.63500229, -0.64072731,... 
                -0.64640252, -0.65203192, -0.65761903, -0.66316698, -0.66867854, -0.67415622, -0.67960224, -0.68501863,... 
                -0.69040721, -0.69576964, -0.70110742, -0.70642194, -0.71171445, -0.71698611, -0.72223799, -0.72747106,...
                -0.73268625, -0.73788439, -0.74306626, -0.74823259, -0.75338406, -0.75852132, -0.76364494, -0.76875551,... 
                -0.77385353, -0.7789395];
        elseif wchoice == 2
            % Ab-initio Data for 3.9 nm fin width (p-type)
            Extra_hole = [0,  0.00020408,  0.00040816,  0.00061224,  0.00081633, 0.00102041,  0.00122449,  0.00142857,...  
                0.00163265,  0.00183673, 0.00204082,  0.0022449 ,  0.00244898,  0.00265306,  0.00285714, 0.00306122,...  
                0.00326531,  0.00346939,  0.00367347,  0.00387755, 0.00408163,  0.00428571,  0.0044898 ,  0.00469388,...  
                0.00489796,  0.00510204,  0.00530612,  0.0055102 ,  0.00571429,  0.00591837, 0.00612245,  0.00632653,...  
                0.00653061,  0.00673469,  0.00693878,  0.00714286,  0.00734694,  0.00755102,  0.0077551 ,  0.00795918,...
                0.00816327,  0.00836735,  0.00857143,  0.00877551,  0.00897959,  0.00918367,  0.00938776,  0.00959184,...  
                0.00979592,  0.01];
            Gate_Bias = [0, -0.44104639, -0.46373664, -0.4789879 , -0.49119328, -0.50172805, -0.51120622, -0.51995458,... 
                -0.52816831, -0.53597334, -0.54345562, -0.5506765 , -0.55768125, -0.56450426, -0.57117227, -0.57770645,... 
                -0.5841239 , -0.59043859, -0.59666208, -0.60280408, -0.60887279, -0.61487519, -0.62081733, -0.62670442,... 
                -0.63254101, -0.63833109, -0.64407819, -0.64978544, -0.65545562, -0.66109121, -0.66669446, -0.67226738,...
                -0.67781181, -0.68332939, -0.68882163, -0.69428992, -0.6997355 , -0.70515953, -0.71056304, -0.71594692,...
                -0.72131203, -0.72665933, -0.73199064, -0.73730504, -0.74260394, -0.74788788, -0.75315743, -0.75841316,... 
                -0.76365558, -0.76888519];
        elseif wchoice == 3
            % Ab-initio Data for 5.4 nm fin width (p-type)
            Extra_hole = [0,  0.00020408,  0.00040816,  0.00061224,  0.00081633, 0.00102041,  0.00122449,  0.00142857,...  
                0.00163265,  0.00183673, 0.00204082,  0.0022449 ,  0.00244898,  0.00265306,  0.00285714, 0.00306122,...  
                0.00326531,  0.00346939,  0.00367347,  0.00387755, 0.00408163,  0.00428571,  0.0044898 ,  0.00469388,...  
                0.00489796,  0.00510204,  0.00530612,  0.0055102 ,  0.00571429,  0.00591837, 0.00612245,  0.00632653,...  
                0.00653061,  0.00673469,  0.00693878,  0.00714286,  0.00734694,  0.00755102,  0.0077551 ,  0.00795918,...
                0.00816327,  0.00836735,  0.00857143,  0.00877551,  0.00897959,  0.00918367,  0.00938776,  0.00959184,...  
                0.00979592,  0.01]; 
            Gate_Bias = [0, -0.4173284 , -0.44055949, -0.45634295, -0.4690758 , -0.48013489, -0.49013216, -0.49939454,... 
                -0.50811706, -0.51642763, -0.52441115, -0.5321289 , -0.53962497, -0.54693602, -0.55408766, -0.56110108,... 
                -0.56799335, -0.57477842, -0.58146785, -0.58807135, -0.5945971 , -0.6010521 , -0.60744236, -0.61377312,... 
                -0.62004892, -0.62627374, -0.63245113, -0.63858421, -0.64467575, -0.65072826, -0.65674397, -0.66272491,... 
                -0.66867288, -0.67458955, -0.6804764 , -0.68633489, -0.69216729, -0.69797131, -0.70375192, -0.70950869,...
                -0.71524216, -0.72095346, -0.72664326, -0.73231269, -0.73796407, -0.74359414, -0.74920576, -0.75479927,... 
                -0.76037521, -0.76593405];
        elseif wchoice == 4
            % Ab-initio Data for 8.7 nm fin width (p-type)
            Extra_hole = [ 0,  0.00020408,  0.00040816,  0.00061224,  0.00081633, 0.00102041,  0.00122449,  0.00142857,...  
                0.00163265,  0.00183673,  0.00204082,  0.0022449 ,  0.00244898,  0.00265306,  0.00285714, 0.00306122,...  
                0.00326531,  0.00346939,  0.00367347,  0.00387755,  0.00408163,  0.00428571,  0.0044898 ,  0.00469388,...  
                0.00489796,  0.00510204,  0.00530612,  0.0055102 ,  0.00571429,  0.00591837,  0.00612245,  0.00632653,...  
                0.00653061,  0.00673469,  0.00693878,  0.00714286,  0.00734694,  0.00755102,  0.0077551 ,  0.00795918,...
                0.00816327,  0.00836735,  0.00857143,  0.00877551,  0.00897959,  0.00918367,  0.00938776,  0.00959184,...  
                0.00979592,  0.01];  
            Gate_Bias = [0, -0.40153703, -0.42562917, -0.44225136, -0.45580189, -0.46765534, -0.47842493, -0.48843797, -0.49789372,... 
                -0.50690743, -0.51557629, -0.52395861, -0.53209998, -0.54003522, -0.54779146, -0.55539046, -0.56284975, -0.57018328,... 
                -0.57740335, -0.58452008, -0.59154212, -0.59847688, -0.60533082, -0.61210958, -0.61881809, -0.6254615 , -0.63204216,... 
                -0.63856533, -0.64503347, -0.6514487 , -0.65781415, -0.66413234, -0.67040546, -0.67663574, -0.68282464, -0.6889741 ,... 
                -0.69508572, -0.70116098, -0.70720127, -0.71320788, -0.71918202, -0.72512461, -0.73103727, -0.73692071, -0.74277592,...
                -0.74860379, -0.75440518, -0.76018092, -0.76593177, -0.77165847];
        elseif wchoice == 5
            % Ab-initio Data for 10.3 nm fin width (p-type)
            
        end
        
    end
    
    % Cleaning up of the Ab-initio Data
    % Removing (0,0) point
    Gate_Bias(1) = [];
    if Simulation_type == 'n'
        Extra_electron(1) = [];
    elseif Simulation_type == 'p'
        Extra_hole(1) = [];
    end
    
    % Procedure to Calculate Surface Bias on Silicon layer based on Displacement Electric Fields being equal at interface
    % Calculation of Displacement Electric Field
    if Simulation_type == 'n'
        Disp_Electric_Field = -(Extra_electron*1.6*10^-19)/(2*3.84*3.84*10^-20);
    elseif Simulation_type == 'p'
        Disp_Electric_Field = -(Extra_hole*1.6*10^-19)/(2*3.84*3.84*10^-20);
    end    
    
    % Calculation of Electric Field in Oxide
    Electric_Field_oxide = Disp_Electric_Field/(3.9 *8.854 * 10^-12);
    % Calculation of Potential Drop across Oxide
    Potential_Drop       = Electric_Field_oxide*10^-9;
    % Calculation of Surface Bias on Silicon layer
    Surface_Bias         = Gate_Bias - Potential_Drop;
    
    % Ab-initio plots
    if Simulation_type == 'n'
        % Plot of Extra Electron vs. Gate Bias
        figure(fileloop+2);
        plot(Gate_Bias,Extra_electron,'o');
        xlabel('Gate Bias(V)','FontSize', 20);
        ylabel('Extra electron(number)','FontSize', 20);
        title('Plot of Extra electron vs. Gate Bias for different methods','FontSize', 20);
        hold on;
    elseif Simulation_type == 'p'
        % Plot of Extra Hole vs. Gate Bias
        figure(fileloop+2);
        plot(Gate_Bias,Extra_hole,'o');
        xlabel('Gate Bias(V)','FontSize', 20);
        ylabel('Extra hole(number)','FontSize', 20);
        title('Plot of Extra hole vs. Gate Bias for different methods','FontSize', 20);
        hold on;
    end


    % Schrodinger-Poisson Solver
    % Choice on whether to run Schrodinger-Poisson solver
    if choice == 'Y'
        % Schrodinger-Poisson Solver uses Gate Bias Directly
        Vg = Gate_Bias;
        % Pre-Initialisation
        if Simulation_type == 'n'
            SP_Extra_electron = ones(1,length(Vg));
        elseif Simulation_type == 'p'
            SP_Extra_hole = ones(1,length(Vg));
        end
        % Loop of Schrodinger-Poisson solvers for Different Gate Biases
        for counter1 = 1:length(Vg)
        % Bias Condition
        OutGenerate.surface_potential =  Vg(counter1);
        if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
            for i = 1:size(OutGenerate.surface.bias,2)
                try
                     % Calling Main File of Effective Mass Schrodinger-Poisson Solver
                     scs(H,OutGenerate,OutGenerate.surface.bias{i}); 

                catch err
                    idSegLast = regexp(err.identifier, '(?<=:)\w+$', 'match');
                    if  strcmp(idSegLast,'ConvergenceFailed')
                        fprintf('\n'); fprintf(err.message); fprintf('\n\n');
                        return
                    else
                        rethrow(err);
                    end
                end
            end

        output_dir  = OutGenerate.control.Output_directory; 
        file_prefix = OutGenerate.control.file_prefix;

        fprintf('\nSaving the results in:\n %s %s \n',output_dir,file_prefix);
        end

        % Processing the Schrodinger Poisson Files
        SP = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '.mat']);
        
        if Simulation_type == 'n'    
            % Using MATLAB Trapezoidal Function to calculate Extra Electron from
            % Schrodinger-Poisson Calculations
            SP_Extra_electron(counter1) = trapz(SP.OutGenerate.x, SP.p_charge - ... 
            SP.n_charge)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);  
        elseif Simulation_type == 'p'
            % Using MATLAB Trapezoidal Function to calculate Extra hole from
            % Schrodinger-Poisson Calculations
            SP_Extra_hole(counter1) = trapz(SP.OutGenerate.x, SP.p_charge - ... 
            SP.n_charge)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);  
        end            
        end  

        if Simulation_type == 'n'
            % Schrodinger-Poisson Plots
            % Plot of Extra Electron vs. Gate Bias
            figure(fileloop+2);
            plot(Gate_Bias,SP_Extra_electron,'--');
            hold on;
            
        elseif Simulation_type == 'p'
            % Schrodinger-Poisson Plots
            % Plot of Extra hole vs. Gate Bias
            figure(fileloop+2);
            plot(Gate_Bias,SP_Extra_hole,'--');
            hold on;
           
        end
    end

    % Density-Gradient Optimization Solvers
    % Processing of Input File for Density-Gradient Solvers
    input_file_DG = files_DG{fileloop};
    OutGenerate=generate(input_file_DG);
    % Density-Gradient Solvers uses Surface Bias 
    Vg = Surface_Bias;
    
    % Optimized Single Fitting Factor Density-Gradient Solver
    for counter5 = 1:length(Vg)
        % Calculating Initial Solution(Poisson) for Density-Gradient Scheme
        OutGenerate.surface_potential =  Vg(counter5);
        [Guess(counter5).Vbetterguess, Guess(counter5).A]   = DG_initial(OutGenerate);
    end
    % Optimization Scheme to find the best-fit Fitting Factor as per least squares method
    fitting_factor_initial_single = 1;
    lb=0;
    ub=100;
    options.TolFun=1e-10;
    options.TolX=1e-25;
    options.MaxIter=500;
    options.MaxFunEvals=2000;
    if Simulation_type == 'n'
        fitting_factor_final_single = lsqnonlin(@abinitio_fitting3,fitting_factor_initial_single,lb,ub,options,OutGenerate,Extra_electron,Vg,Guess);
    elseif Simulation_type == 'p'
        fitting_factor_final_single = lsqnonlin(@abinitio_fitting4,fitting_factor_initial_single,lb,ub,options,OutGenerate,Extra_hole,Vg,Guess);
    end
    fprintf('DG fitting parameter = %f',fitting_factor_final_single); 
    
    % Pre-Initialisation
    if Simulation_type == 'n'
        DG_Extra_electron_single = ones(1,length(Vg));
    elseif Simulation_type == 'p'
        DG_Extra_hole_single = ones(1,length(Vg));
    end
    % Calculaion of Extra Electron as per Optimized Fitting Factor
    for counter = 1:length(Vg)
        OutGenerate.surface_potential =  Vg(counter);
        [Vbetterguess, A]   = DG_initial(OutGenerate);
        [~] = scsdg(fitting_factor_final_single,OutGenerate,Vbetterguess,A);
        % Processing the Density-Gradient files
        DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter)) '_DG.mat']);
        if Simulation_type == 'n'
            % Using MATLAB Trapezoidal Function to calculate Extra Electron from
            % Density-Gradient Calculations
            DG_Extra_electron_single(counter) = trapz(DG.x, DG.p_charge_DG - ... 
            DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);
        elseif Simulation_type == 'p'
            % Using MATLAB Trapezoidal Function to calculate Extra hole from
            % Density-Gradient Calculations
            DG_Extra_hole_single(counter) = trapz(DG.x, DG.p_charge_DG - ... 
            DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);
        end
    end
    
    if Simulation_type == 'n'
        % Density-Gradient Single Fitting Factor plots 
        % Plot of Extra Electron vs. Gate Bias
        figure(fileloop+2);
        plot(Gate_Bias,DG_Extra_electron_single,'-.');
        hold on;
    elseif Simulation_type == 'p'
        % Density-Gradient Single Fitting Factor plots 
        % Plot of Extra hole vs. Gate Bias
        figure(fileloop+2);
        plot(Gate_Bias,DG_Extra_hole_single,'-.');
        hold on;
    end
    
    % Plot of Fitting Factor vs. Applied Gate Bias
    figure(1);
    if Simulation_type == 'n'
        hline = refline(0,fitting_factor_final_single);
    elseif Simulation_type == 'p'
        hline = refline(0,fitting_factor_final_single);
    end
    % Colour Scheme of Plot
    if fileloop == 1
        hline.Color = 'b';
    elseif fileloop == 2
        hline.Color = 'g';
    elseif fileloop == 3
        hline.Color = 'r';
    elseif fileloop == 4
        hline.Color = 'c';
    elseif fileloop == 5
        hline.Color = 'm';
    end
    hline.LineStyle = '--';
    xlabel('Gate Bias(V)','FontSize', 20);
    ylabel('DG Fitting Factor','FontSize', 20);
    title('Plot of Variation of Density-Gradient Fitting Factor with Applied Gate Voltage','FontSize', 20);
    hold on;

    % Optimized Multiple Fitting Factor Density-Gradient Solver
    if strcmp(Simulation_procedure,'optimization')    
        % Pre-Initialisation
        if Simulation_type == 'n'
            DG_Extra_electron = ones(1,length(Vg));
        elseif Simulation_type == 'p'
            DG_Extra_hole = ones(1,length(Vg));
        end
        fitting_factor_final = ones(1,length(Vg));
        % Loop of Density-Gradient solver for finding the optimized fitting factor at that Gate Bias
        for counter1 = 1:length(Vg)   
            % Calculating Initial Solution(Poisson) for Density-Gradient Scheme
            OutGenerate.surface_potential =  Vg(counter1);
            [Vbetterguess, A]   = DG_initial(OutGenerate);
            
            % Optimization Scheme to find the best Fitting Factor at current Gate Bias
            fitting_factor_initial = fitting_factor_final_single;
            lb=0;
            ub=100;
            options.TolFun=1e-10;
            options.TolX=1e-15;
            options.MaxIter=500;
            options.MaxFunEvals=2000;
            if Simulation_type == 'n'
                fitting_factor_final(counter1) = lsqnonlin(@abinitio_fitting2,fitting_factor_initial,lb,ub,options,OutGenerate,Extra_electron(counter1),Vbetterguess,A);
            elseif Simulation_type == 'p'
                fitting_factor_final(counter1) = lsqnonlin(@abinitio_fitting5,fitting_factor_initial,lb,ub,options,OutGenerate,Extra_hole(counter1),Vbetterguess,A);
            end
            fprintf('DG fitting parameter = %f',fitting_factor_final); 

            % Processing Density-Gradient files
            DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '_DG.mat']);
            if Simulation_type == 'n'
                % Using MATLAB Trapezoidal Function to calculate Extra Electron from
                % Density-Gradient Calculations
                DG_Extra_electron(counter1) = trapz(DG.x, DG.p_charge_DG - ... 
                DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);
            elseif Simulation_type == 'p'
                % Using MATLAB Trapezoidal Function to calculate Extra hole from
                % Density-Gradient Calculations
                DG_Extra_hole(counter1) = trapz(DG.x, DG.p_charge_DG - ... 
                DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);
            end
        
        end
        
        % Plots of Multiple Fitting Factor Density-Gradient Solver
        if Simulation_type == 'n'
            % Plot of Extra electron vs. Gate Bias
            figure(fileloop+2);
            plot(Gate_Bias,DG_Extra_electron,'-','LineWidth',2);
            hold on;
        elseif Simulation_type == 'p'
            % Plot of Extra hole vs. Gate Bias
            figure(fileloop+2);
            plot(Gate_Bias,DG_Extra_hole,'-','LineWidth',2);
            hold on;
        end
        % Figure Specifications
        figure(fileloop+2);
        if choice == 'Y'
            legend('Ab-initio','Effective Mass Schrodinger-Poisson','Density-Gradient(Single Fit factor)','Density-Gradient(Multiple Fit-factor)');
        else
            legend('Ab-initio','Density-Gradient(Single Fit factor)','Density-Gradient(Multiple Fit-factor)');
        end

        % Plotting Variation of Fitting Factor vs. Applied Gate Bias
        figure(1);        
        if fileloop == 1
            p1 = plot(Gate_Bias,fitting_factor_final,'bx-','LineWidth',2 ); 
        elseif fileloop == 2
            p2 = plot(Gate_Bias,fitting_factor_final,'gx-','LineWidth',2 ); 
        elseif fileloop == 3
            p3 = plot(Gate_Bias,fitting_factor_final,'rx-','LineWidth',2 ); 
        elseif fileloop == 4
            p4 = plot(Gate_Bias,fitting_factor_final,'cx-','LineWidth',2 ); 
            p4d = plot(Gate_Bias,fitting_factor_final,'cx-','LineWidth',2 ); 
        elseif fileloop == 5
            p5 = plot(Gate_Bias,fitting_factor_final,'mx-','LineWidth',2 );
            p5d = plot(Gate_Bias,fitting_factor_final,'mx-','LineWidth',2 );
        end       

        % Error Calculations  
        if Simulation_type == 'n'            
            % Plot of Error between Single and Multiple Fitting Factors and Ab-initio
            Percent_Error_single = abs(((DG_Extra_electron_single - Extra_electron)./Extra_electron)*100);
            figure(8)
            if fileloop == 1
                t1 = plot(Gate_Bias,Percent_Error_single,'b--','LineWidth',2);
            elseif fileloop == 2
                t3 = plot(Gate_Bias,Percent_Error_single,'g--','LineWidth',2);
            elseif fileloop == 3
                t5 = plot(Gate_Bias,Percent_Error_single,'r--','LineWidth',2); 
            elseif fileloop == 4
                t7 = plot(Gate_Bias,Percent_Error_single,'c--','LineWidth',2);
            elseif fileloop == 5
                t9 = plot(Gate_Bias,Percent_Error_single,'m--','LineWidth',2); 
            end
            xlabel('Gate Bias(V)','FontSize', 20);
            ylabel('Percent Error','FontSize', 20);
            title('Percent Error between Single and Multiple Fitting Factos in Density-Gradient','FontSize', 20);
            hold on;
        elseif Simulation_type == 'p'
            % Plot of Error between Single and Multiple Fitting Factors and Ab-initio
            Percent_Error_single = abs(((DG_Extra_hole_single - Extra_hole)./Extra_hole)*100);
            figure(8)
            if fileloop == 1
                t1 = plot(Gate_Bias,Percent_Error_single,'b--','LineWidth',2);
            elseif fileloop == 2
                t3 = plot(Gate_Bias,Percent_Error_single,'g--','LineWidth',2);
            elseif fileloop == 3
                t5 = plot(Gate_Bias,Percent_Error_single,'r--','LineWidth',2); 
            elseif fileloop == 4
                t7 = plot(Gate_Bias,Percent_Error_single,'c--','LineWidth',2);
            elseif fileloop == 5
                t9 = plot(Gate_Bias,Percent_Error_single,'m--','LineWidth',2); 
            end
            xlabel('Gate Bias(V)','FontSize', 20);
            ylabel('Percent Error','FontSize', 20);
            title('Percent Error between Single and Multiple Fitting Factos in Density-Gradient','FontSize', 20);
            hold on;
        end

    elseif strcmp(Simulation_procedure,'sweep')
        % Defining the sweep range of fitting factors
        fitting_factor_final = linspace(0 , 0.9 , 10 );
        for counter2 = 1:length(fitting_factor_final)
            for counter1 = 1:length(Gate_Bias)
                % Calculating Initial Solution(Poisson) for Density-Gradient Scheme
                OutGenerate.surface_potential =  Vg(counter1);
                [Vbetterguess, A]   = DG_initial(OutGenerate);              
                [~] = scsdg(fitting_factor_final(counter2),OutGenerate,Vbetterguess,A);

                % Processing the Density Gradient Files
                DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '_DG.mat']);
                
                if Simulation_type == 'n'
                    % Using MATLAB Trapezoidal Function to calculate Extra Electron from
                    % Density-Gradient Calculations
                    DG_Extra_electron(counter1) = trapz(DG.x, DG.p_charge_DG - ... 
                    DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16); 
                elseif Simulation_type == 'p'
                    % Using MATLAB Trapezoidal Function to calculate Extra hole from
                    % Density-Gradient Calculations
                    DG_Extra_hole(counter1) = trapz(DG.x, DG.p_charge_DG - ... 
                    DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16); 
                end
            end     
            
        % Plots of Sweep Fitting Factors
        if Simulation_type == 'n'
            % Plot of Extra electron vs. Gate Bias
            figure(fileloop+2);
            plot(Gate_Bias,DG_Extra_electron,'-','LineWidth',2);
            hold on;    
        elseif Simulation_type == 'p'
            % Plot of Extra hole vs. Gate Bias
            figure(fileloop+2);
            plot(Gate_Bias,DG_Extra_hole,'-','LineWidth',2);
            hold on;  
        end  
        end

    % Plot Specifications  
    if choice == 'Y'
        figure(fileloop+2);        
        legend('ab-initio','Effective Mass Schrodinger-Poisson','Density-Gradient (Single Fit factor)','fitting factor = 0',...
        'fitting factor = 1', 'fitting factor = 2', 'fitting factor = 3', 'fitting factor = 4', 'fitting factor = 5',... 
        'fitting factor = 6', 'fitting factor = 7', 'fitting factor = 8' ,'fitting factor = 9' ,'Location','SouthWest');
    else
        figure(fileloop+2)
        legend('ab-initio','Density-Gradient (Single Fit factor)','fitting factor = 0',...
        'fitting factor = 1', 'fitting factor = 2', 'fitting factor = 3', 'fitting factor = 4', 'fitting factor = 5',... 
        'fitting factor = 6', 'fitting factor = 7', 'fitting factor = 8' ,'fitting factor = 9' ,'Location','NorthEast');
    end
    end    
end

if strcmp(Simulation_procedure,'optimization')
    % Figure 1 Specifications
    figure(1)
    legend([p1 p2 p3 p4 hline p4d],'3.2nm','3.9nm','5.4nm','8.7nm','Single Fitting Factor','Multiple Fitting Factor');
    hold on;

    % Figure 8 Specifications
    figure(8)
    legend([t1 t3 t5 t7],{'3.2nm','3.9nm','5.4nm','8.7nm','10.3nm'});
    hold on;
end
end

