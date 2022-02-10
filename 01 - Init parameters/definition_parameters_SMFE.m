%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal cutting condition selection for high quality receptance         
% measurements by Sweep Milling Force Excitation
% 
% Authors: Oier Franco, Xavier Beudaert, Alex Iglesias, Zoltan Dombovari,
%          Kaan Erkorkmaz, Jokin Munoa. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definition of the tool
tool.Z      = 2;            % number of teeth of the tool
tool.D      = 80/1000;      % diameter of the tool (m)
tool.kappa  = 45*pi/180;    % angle kappa of the insert (rad)

% Variable pitch
tool.VariablePitch.Enable='false';
tool.VariablePitch.Angles=[0 44 90 134 180 224 270 314 360]*pi/180; % angular position of the teeth (rad)

% run out
tool.RunOut.Enable='false';
tool.RunOut=zeros(tool.Z,3);

%% Defition of the process
process.fz   = 0.1/1000;     % feed per tooth (m/tooth)
process.ap   = 2/1000;       % axial depth of cut (m)
process.ae   = 12/1000;      % width of cut or radial engagement (m)
process.Nrpm = 600;          % spindle speed (rpm)
process.N    = process.Nrpm*2*pi/60; % spindle speed (rad per second)

% cutting coefficients

process.cutting_coef.Kt  = 1889.2*1e6; % tangential cutting coefficient (N/m2)
process.cutting_coef.Kr  = 806.1*1e6 ; % radial cutting coefficient (N/m2)
process.cutting_coef.Ka  = 291.0*1e6;  % axial cutting coefficient (N/m2)

process.cutting_coef.Kte = 63.12*1e3;  % tangential edge coefficient (N/m)
process.cutting_coef.Kre = 113.7*1e3;  % radial edge coefficient (N/m)
process.cutting_coef.Kae = -3.1*1e3;   % axial edge coefficient (N/m)

% Sweep Milling Force Excitation
process.SMFE.Enable       = 'true';
process.SMFE.Acceleration = 1.5;       % Acceleration rate for spindle speed modification
% Spindle speed variation
process.SSV.Enable        = 'false';
process.SSV.Amplitude     = 10;        % Oscillation amplitude (% of nominal spindle speed)
process.SSV.Frequency     = 1;         % Frequency of the sine oscillation (Hz)


% calulation of the angles of the cut
if strcmp(process.operation_type,'up_milling')
    process.start_angle=0; % angle phi from which the tooth is in cut (rad)
    process.exit_angle=acos((0.5*tool.D - process.ae)/(0.5*tool.D)); % angle phi after which the tooth is out of cut (rad)
elseif strcmp(process.operation_type,'down_milling')
    process.start_angle=acos((process.ae - 0.5*tool.D)/(0.5*tool.D));
    process.exit_angle=pi;
elseif strcmp(process.operation_type,'slotting')
    process.start_angle=acos((0.5*process.ae) / (0.5*tool.D));
    process.exit_angle=pi-process.start_angle;
elseif strcmp(process.operation_type,'full_immersion')
    process.start_angle=0;
    process.exit_angle=pi;
else
    error('Unknown process.operation_type')
end


%% Definition of the machine
% Modal parameters
machine.dynamic.number_of_mode=3; % number of modes
machine.dynamic.fn=[33.37, 57.44, 58.61]; % natural frequencies (Hz)
machine.dynamic.wn=2*pi*machine.dynamic.fn; % natural frequencies (rad/s)
machine.dynamic.K=[45.5, 32.1, 36.7]*1e6; % rigidez (N/m)
machine.dynamic.mass=machine.dynamic.K ./ (machine.dynamic.wn.^2); % reflected mass (kg)
machine.dynamic.ksi=[0.061 0.0822 0.0287]; % damping ratio ()
machine.dynamic.c=2*machine.dynamic.mass .* machine.dynamic.ksi .* machine.dynamic.wn; % modal damping (Ns/m)
machine.dynamic.modal_vectors=[-0.055 0.0199 0.999 ;
                               +0.6245 0.9968 4.0028e-4;
                               0.7724 -0.0161 0.0096 ];

% Rotation from the Machine Coordinate System to the Tool Coordinate System                
machine.angles.alpha=0; 
machine.angles.psi=0;
machine.angles.theta=0;
                
% Rotation from the Machine Coordinate System to the Tool Coordinate System
machine.ROT_MCS_TCS = [ cos(machine.angles.alpha)*cos(machine.angles.psi)-sin(machine.angles.alpha)*cos(machine.angles.theta)*sin(machine.angles.psi),  cos(machine.angles.alpha)*sin(machine.angles.psi)+sin(machine.angles.alpha)*cos(machine.angles.theta)*cos(machine.angles.psi), sin(machine.angles.alpha)*sin(machine.angles.theta);...
                       -sin(machine.angles.alpha)*cos(machine.angles.psi)-cos(machine.angles.alpha)*cos(machine.angles.theta)*sin(machine.angles.psi), -sin(machine.angles.alpha)*sin(machine.angles.psi)+cos(machine.angles.alpha)*cos(machine.angles.theta)*cos(machine.angles.psi), cos(machine.angles.alpha)*sin(machine.angles.theta);...
                        sin(machine.angles.theta)*sin(machine.angles.psi),                                                                             -sin(machine.angles.theta)*cos(machine.angles.psi),                                                                                          cos(machine.angles.theta)];                        
% Rotation from the Tool Coordinate System to the Machine Coordinate System
machine.ROT_TCS_MCS=inv(machine.ROT_MCS_TCS);

% Frequency discretization of the FRF
machine.frequency_vector=linspace(0,200,1000).'; % frequency vector (Hz)
machine.FRF=zeros(size(machine.dynamic.modal_vectors,1),size(machine.dynamic.modal_vectors,1),length(machine.frequency_vector));
for IdxMode=1:machine.dynamic.number_of_mode % loop over modes
    for IdxFRF1=1:size(machine.FRF,1)
        for IdxFRF2=1:size(machine.FRF,2)
            machine.FRF(IdxFRF1,IdxFRF2,:)= machine.FRF(IdxFRF1,IdxFRF2,:)+ ...
                machine.dynamic.modal_vectors(IdxFRF1,IdxMode)*machine.dynamic.modal_vectors(IdxFRF2,IdxMode) * ...
                permute(1/machine.dynamic.mass(IdxMode) ./((2*pi*machine.dynamic.fn(IdxMode))^2 - (2*pi*machine.frequency_vector).^2 ...
                                         +2*1i*machine.dynamic.ksi(IdxMode)*(2*pi*machine.dynamic.fn(IdxMode))*(2*pi*machine.frequency_vector))  ,[3 2 1]);
        end
    end
end
%  Plot FRF
figure; hold on; box on; grid on
for IdxFRF1=1:size(machine.FRF,1)
    for IdxFRF2=1:size(machine.FRF,2)
        plot(machine.frequency_vector,abs(squeeze(machine.FRF(IdxFRF1,IdxFRF2,:))),'DisplayName',strcat('FRF',num2str(IdxFRF1),num2str(IdxFRF2)));
    end
end
xlabel('Frequency (Hz)'); ylabel('amplitude (m / N)'); title('Defined machine tool dynamics')
legend('show')
process.ROT_MCS_TCS=machine.ROT_MCS_TCS;

%% Definition of the workpiece
workpiece.dynamic.number_of_mode=0; % number of modes
workpiece.dynamic.fn=[50 30 20]; % natural frequencies (Hz)
workpiece.dynamic.wn=2*pi*workpiece.dynamic.fn; % natural frequencies (rad/s)
workpiece.dynamic.ksi=[5 3 4 1]/100; % damping
workpiece.dynamic.K=[28 20 30]*1e6; % rigidez (N/m)
workpiece.dynamic.modal_vectors=[ 0.1822, 0.2094, 0.7853;...
                               -0.5740, 0.7507, 0.6163;...
                               -0.7983, 0.6266, -0.0598];

% Frequency discretization of the FRF
workpiece.frequency_vector=linspace(0,200,200).'; % frequency vector (Hz)
workpiece.FRF=zeros(size(workpiece.dynamic.modal_vectors,1),size(workpiece.dynamic.modal_vectors,1),length(workpiece.frequency_vector));
for IdxMode=1:workpiece.dynamic.number_of_mode % loop over modes
    for IdxFRF1=1:size(workpiece.FRF,1)
        for IdxFRF2=1:size(workpiece.FRF,2)
            workpiece.FRF(IdxFRF1,IdxFRF2,:)= workpiece.FRF(IdxFRF1,IdxFRF2,:)+ ...
                workpiece.dynamic.modal_vectors(IdxFRF1,IdxMode)*workpiece.dynamic.modal_vectors(IdxFRF2,IdxMode) * ...
                permute(1/workpiece.dynamic.mass(IdxMode) ./((2*pi*workpiece.dynamic.fn(IdxMode))^2 - (2*pi*workpiece.frequency_vector).^2 ...
                                         +2*1i*workpiece.dynamic.ksi(IdxMode)*(2*pi*workpiece.dynamic.fn(IdxMode))*(2*pi*workpiece.frequency_vector))  ,[3 2 1]);
        end
    end
end

%% Definition of the solver
solver.number_of_tool_revolutions=650;  % Acceleration
solver.enable_out_of_cut_non_linearity=1;
solver.check_cut_of_N_previous_teeth=0; % tool.Z, check if the Nth previous teeth were cutting more material

solver.max_dt=0.0005; % solver maximum time step (s)
solver.max_angular_step=1*pi/180; % solver maximum angular step
solver.min_number_points_of_highest_frequency=20; % minimum number of sample for the highest frequency in the model

solver.method='Euler';

% adjustement of solver parameters depending on the spindle speed N
temp_dt=min([solver.max_dt,...
    1/(solver.min_number_points_of_highest_frequency * max([machine.dynamic.fn workpiece.dynamic.fn])), ... % condition for the highest modal frequency
    solver.max_angular_step/process.N]);

solver.number_of_angular_positions_between_2_teeth=ceil( ( (2*pi/tool.Z)) /(process.N*temp_dt)); % how many simulation steps between 2 teeth
solver.angular_step= (2*pi/tool.Z) /solver.number_of_angular_positions_between_2_teeth; % angular discretization (rad)
solver.dt=solver.angular_step/process.N; % time discretization (s)
solver.time=(0:solver.dt:(2*pi/process.N)*solver.number_of_tool_revolutions)';
solver.number_of_simulation_steps=length(solver.time);

clear temp_dt
