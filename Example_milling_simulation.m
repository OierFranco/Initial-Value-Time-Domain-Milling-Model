%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal cutting condition selection for high quality receptance         
% measurements by Sweep Milling Force Excitation
% 
% Authors: Oier Franco, Xavier Beudaert, Alex Iglesias, Zoltan Dombovari,
%          Kaan Erkorkmaz, Jokin Munoa. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear variables
addpath('00 - Functions');
addpath('01 - Init parameters');

%%
% Define the type of cutting operation 
process.operation_type='down_milling'; % Options:    down_milling
                                       %             up_milling 
                                       %             slotting 

% Load pre-defined parameters:

%     definition_parameters_SMFE;
    definition_parameters_regular_cut;
    
[tool,workpiece]=milling_time_domain_simulation(tool,process,machine,workpiece,solver);

%% Save simulation data for next step:
% Format: X Y Z accel + X Y Z forces + Time
    Force = tool.forces.cutting_forces_MCS;
    Accel = tool.vib_acc_MCS;
    Time = solver.time;

if strcmp(process.operation_type,'up_milling')
    Test_up = [Accel(:,1), Accel(:,2), Accel(:,3), Force(:,1), Force(:,2), Force(:,3),Time];
elseif   strcmp(process.operation_type,'down_milling') 
    Test_down = [Accel(:,1), Accel(:,2), Accel(:,3), Force(:,1), Force(:,2), Force(:,3),Time];
else 
    Test_central = [Accel(:,1), Accel(:,2), Accel(:,3), Force(:,1), Force(:,2), Force(:,3),Time];
end 

% From Workspace, Save as .mat, Test_up, Test_down or Test_central results.