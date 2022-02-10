%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal cutting condition selection for high quality receptance         
% measurements by Sweep Milling Force Excitation
% 
% Authors: Oier Franco, Xavier Beudaert, Alex Iglesias, Zoltan Dombovari,
%          Kaan Erkorkmaz, Jokin Munoa. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   The inputs are defined in the file "definition_parameters.m"
%
% OUTPUT:
%   - Cutting forces
%   - Tool-Workpiece vibration (optional)

function [tool,workpiece]=milling_time_domain_simulation(tool,process,machine,workpiece,solver)

%% Pre allocation
% angles
if strcmp(process.SMFE.Enable,'true') % SMFE process activated
    
    tool_angle = process.N * solver.time + 0.5*process.SMFE.Acceleration .* solver.time.*solver.time;
    Velocity_rad_s = process.N + process.SMFE.Acceleration.*solver.time;
    Velocity_rpm   = Velocity_rad_s.*(60)./(2*pi);

    figure(2)
    subplot(2,1,1)
    plot(solver.time,Velocity_rpm,'Linewidth',2,'DisplayName', 'SMFE');
    hold on
    plot([solver.time(1), solver.time(end)],[Velocity_rpm(1),Velocity_rpm(1)],'k','Linewidth',2,'DisplayName', 'Regular');
    xlim([solver.time(1), solver.time(end)]);
    ylim([0 max(Velocity_rpm)*1.5])
    xlabel('Time (s)'); ylabel('Spindle speed (rpm)')
    grid on


elseif strcmp(process.SSV.Enable,'true') % Spindle Speed Variation activated
    
    SSV_sine_wave   = process.N * (process.SSV.Amplitude/100) * sin(2*pi*process.SSV.Frequency * solver.time);
    tool_angle      = process.N * solver.time + SSV_sine_wave; % angular position of the tool

else % constant spindle speed
    
    tool_angle = process.N * solver.time;

end

tooth_psi_angle=zeros(solver.number_of_simulation_steps,tool.Z); % angular position of the tooth [0 2*pi]
for IdxTooth=1:tool.Z
    if strcmp(tool.VariablePitch.Enable,'true') % Variable pitch activated
        tooth_psi_angle(:,IdxTooth)=wrapTo2Pi(tool_angle(:)-tool.VariablePitch.Angles(IdxTooth)); % angular position of the tooth [0 2*pi]
    else
        tooth_psi_angle(:,IdxTooth)=wrapTo2Pi(tool_angle(:)-(IdxTooth-1) * (2*pi/tool.Z) ); % angular position of the tooth [0 2*pi]
    end
end

% cutting area
h_total=zeros(solver.number_of_simulation_steps,tool.Z); % total chip tickness (static + dynamic)
depth_of_cut=zeros(solver.number_of_simulation_steps,1);

% Forces
tool.forces.cutting_forces_teeth_TCS=zeros(solver.number_of_simulation_steps,tool.Z,3); % force of each tooth in the TCS
tool.forces.cutting_forces_TCS=zeros(solver.number_of_simulation_steps,3);
tool.forces.cutting_forces_MCS=zeros(solver.number_of_simulation_steps,3);
cutting_forces_modal_space=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode); % Modal Coordinates dim: 1

% tool center complete motion (initial movement + vibration )in the Tool Coordinate System
tool.pos_TCS=zeros(solver.number_of_simulation_steps,3); % position of the tool center in the TCS (m)
tool.vel_TCS=zeros(solver.number_of_simulation_steps,3); % velocity of the tool center in the TCS (m/s)
tool.acc_TCS=zeros(solver.number_of_simulation_steps,3); % acceleration of the tool center in the TCS (m/s^2)

% tool center vibration in the Machine Coordinate System
tool.vib_pos_MCS=zeros(solver.number_of_simulation_steps,3); % vibration of the tool center in the MCS (m)
tool.vib_vel_MCS=zeros(solver.number_of_simulation_steps,3);
tool.vib_acc_MCS=zeros(solver.number_of_simulation_steps,3);

% tool center movement in the Machine Coordinate System
tool.pos_MCS=zeros(solver.number_of_simulation_steps,3); % position of the tool center in the MCS (m)
tool.vel_MCS=zeros(solver.number_of_simulation_steps,3);
tool.acc_MCS=zeros(solver.number_of_simulation_steps,3);

% workpiece vibration in the Machine Coordinate System
workpiece.vib_pos_MCS=zeros(solver.number_of_simulation_steps,3);
workpiece.vib_vel_MCS=zeros(solver.number_of_simulation_steps,3);
workpiece.vib_acc_MCS=zeros(solver.number_of_simulation_steps,3);

% tool center movement in the Machine Coordinate System
tool.mvt_pos_MCS=zeros(solver.number_of_simulation_steps,3); % movement of the tool center in the MCS (m)
tool.mvt_vel_MCS=zeros(solver.number_of_simulation_steps,3);
tool.mvt_acc_MCS=zeros(solver.number_of_simulation_steps,3);

% Tool movement without vibration
tool_start_pos=[-tool.D/2 - 0.005, 0, -process.ap];

tool.mvt_pos_MCS(:,1) = tool_start_pos(1) + process.fz * tool.Z * (process.N/(2*pi)) * solver.time; % constant feed per tooth without control;
tool.mvt_pos_MCS(:,2) = tool_start_pos(2);
tool.mvt_pos_MCS(:,3) = tool_start_pos(3);

% Initilization first index of tool.pos_TCS
tool.pos_TCS(1,:) =  tool.mvt_pos_MCS(1,:);

% modal motion in the Modal Space
acc_modes_modal=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode);
vel_modes_modal=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode);
pos_modes_modal=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode);

% modal motion in the Machine Coordinate System
acc_modes_MCS=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode,3);
vel_modes_MCS=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode,3);
pos_modes_MCS=zeros(solver.number_of_simulation_steps,machine.dynamic.number_of_mode,3);

second_tooth_has_cut='false';

%% Main loop with fixed time step
for IdxStep=1:solver.number_of_simulation_steps-1 % loop over all the tool positions
    %% compute the cutting forces
    for IdxTooth=1:tool.Z % loop over all the teeth
        if tooth_psi_angle(IdxStep,IdxTooth)>=process.start_angle && tooth_psi_angle(IdxStep,IdxTooth)<=process.exit_angle
            
            if IdxTooth==2 && IdxStep~=1
                second_tooth_has_cut='true';
            end
            % Compute the total chip thickness
            if strcmp(second_tooth_has_cut,'false') % used only at the very beginning
                % it is the first tooth cutting so the nominal chip tickness is used
                h_total(IdxStep,IdxTooth) = process.fz * sin(tooth_psi_angle(IdxStep,IdxTooth)) * sin(tool.kappa);
            else
                % position difference between the actual and the previous tool positions at this tooth angle (in the MCS)
                previous_tooth=mod(IdxTooth-1-1,tool.Z)+1;
                magic_percentage=0.2*1; % ref 0.2 to look around the index at which the previous tooth was cutting
                Idx_to_look=max( 1, IdxStep-floor(solver.number_of_angular_positions_between_2_teeth*(1+magic_percentage))) ...
                    : min(solver.number_of_simulation_steps-1 , IdxStep-ceil(solver.number_of_angular_positions_between_2_teeth*(1-magic_percentage)));
                % Interpolation of 2*pi periodic function
                unwrap_angles = unwrap([tooth_psi_angle(Idx_to_look,previous_tooth)',tooth_psi_angle(IdxStep,IdxTooth)]);
                previous_tool_pos_at_this_angle = interp1(...
                    unwrap_angles(1:end-1) ...
                    ,tool.pos_TCS(Idx_to_look,:) ...
                    ,unwrap_angles(end));              
                real_tool_delta = (tool.pos_TCS(IdxStep,:) - tool.pos_TCS(IdxStep-solver.number_of_angular_positions_between_2_teeth,:))';
                % position difference between the actual and the previous workpiece positions at this tooth angle (in the MCS)
                real_workpiece_delta = (workpiece.vib_pos_MCS(IdxStep,:) - workpiece.vib_pos_MCS(IdxStep-solver.number_of_angular_positions_between_2_teeth,:))';
                % relative displacement between the tool and workpiece (in the MCS)
                delta_r = real_tool_delta - real_workpiece_delta;
                
                % projection in the direction of the chip thickness
                u = [sin(tooth_psi_angle(IdxStep,IdxTooth))*sin(tool.kappa); cos(tooth_psi_angle(IdxStep,IdxTooth))*sin(tool.kappa); -cos(tool.kappa)];
                % projection of the radial runout on the chip thickness direction
                u_runout = [sin(tooth_psi_angle(IdxStep,IdxTooth)) ; cos(tooth_psi_angle(IdxStep,IdxTooth)) ; 0];
                h_total(IdxStep,IdxTooth) = u' * delta_r + u_runout' * tool.RunOut(IdxTooth,:)';
            end
            
            depth_of_cut(IdxStep) = workpiece.vib_pos_MCS(IdxStep,3) - tool.pos_TCS(IdxStep,3); % NOT VIB but workpice geometry ???
            
            % out of cut nonlinearity if the chip thickness is negative
            if h_total(IdxStep,IdxTooth)<0
                if solver.enable_out_of_cut_non_linearity==1
                    h_total(IdxStep,IdxTooth) = 0;
                end
            end
            
            % compute the cutting forces of the tooth in local coordinates
            % with linear cutting force model + edge effect
            ft = (process.cutting_coef.Kt * depth_of_cut(IdxStep) / sin(tool.kappa)) * h_total(IdxStep,IdxTooth) ...
                + process.cutting_coef.Kte* depth_of_cut(IdxStep) / sin(tool.kappa)  ;
            fr = (process.cutting_coef.Kr * depth_of_cut(IdxStep) / sin(tool.kappa)) * h_total(IdxStep,IdxTooth) ...
                + process.cutting_coef.Kre* depth_of_cut(IdxStep) / sin(tool.kappa)  ;
            fa = (process.cutting_coef.Ka * depth_of_cut(IdxStep) / sin(tool.kappa)) * h_total(IdxStep,IdxTooth) ...
                + process.cutting_coef.Kae* depth_of_cut(IdxStep) / sin(tool.kappa)  ;
            
            % projection of the cutting force of the tooth in the Tool Coordinates System
            ROT_LCS_TCS = [-cos(tooth_psi_angle(IdxStep,IdxTooth))    -sin(tool.kappa)*sin(tooth_psi_angle(IdxStep,IdxTooth))    -cos(tool.kappa)*sin(tooth_psi_angle(IdxStep,IdxTooth));
                            sin(tooth_psi_angle(IdxStep,IdxTooth))    -sin(tool.kappa)*cos(tooth_psi_angle(IdxStep,IdxTooth))    -cos(tool.kappa)*cos(tooth_psi_angle(IdxStep,IdxTooth));
                            0                                          cos(tool.kappa)                                           -sin(tool.kappa)];           
            tool.forces.cutting_forces_teeth_TCS(IdxStep,IdxTooth,:) = ROT_LCS_TCS * [ft;fr;fa];
        end
    end % END loop over the teeth
    
    % sum the forces given by all the teeth
    tool.forces.cutting_forces_TCS(IdxStep,1:3) = sum(tool.forces.cutting_forces_teeth_TCS(IdxStep,:,:),2);
    
    % projection of the total cutting forces in the Machine Coordinates System
    tool.forces.cutting_forces_MCS(IdxStep,1:3) = machine.ROT_TCS_MCS * tool.forces.cutting_forces_TCS(IdxStep,1:3)';
    
    %% solve the differential equation for the tool side 
    for mode=1:machine.dynamic.number_of_mode
        % transformation of the force to modal space
        cutting_forces_modal_space(IdxStep,mode) = machine.dynamic.modal_vectors(:,mode)' * tool.forces.cutting_forces_MCS(IdxStep,1:3)';
        
        % modal dynamical equation
        acc_modes_modal(IdxStep+1,mode) = (cutting_forces_modal_space(IdxStep,mode)...
                                     - machine.dynamic.c(mode) * vel_modes_modal(IdxStep,mode)...
                                     - machine.dynamic.K(mode) * pos_modes_modal(IdxStep,mode)) ...
                                    / machine.dynamic.mass(mode) ;
                
        if strcmp(solver.method,'Euler')
            vel_modes_modal(IdxStep+1,mode) = vel_modes_modal(IdxStep,mode) + acc_modes_modal(IdxStep+1,mode) * solver.dt;
            pos_modes_modal(IdxStep+1,mode) = pos_modes_modal(IdxStep,mode) + vel_modes_modal(IdxStep+1,mode) * solver.dt;
        else
            error('unkown solver.method')
        end
        
        % Projection to the Machine Coordinate System
        acc_modes_MCS(IdxStep+1,mode,:) = machine.dynamic.modal_vectors(:,mode) * acc_modes_modal(IdxStep+1,mode);
        vel_modes_MCS(IdxStep+1,mode,:) = machine.dynamic.modal_vectors(:,mode) * vel_modes_modal(IdxStep+1,mode);
        pos_modes_MCS(IdxStep+1,mode,:) = machine.dynamic.modal_vectors(:,mode) * pos_modes_modal(IdxStep+1,mode);
    end % END loop over the modes
    
    % Sum the effect of all the modes in the Machine Coordinate System
    tool.vib_acc_MCS(IdxStep+1,:) = sum(acc_modes_MCS(IdxStep+1,:,:),2);
    tool.vib_vel_MCS(IdxStep+1,:) = sum(vel_modes_MCS(IdxStep+1,:,:),2);
    tool.vib_pos_MCS(IdxStep+1,:) = sum(pos_modes_MCS(IdxStep+1,:,:),2);
    
    % Superimpose the tool movement to the tool vibration
    tool.acc_MCS(IdxStep+1,:) = tool.vib_acc_MCS(IdxStep+1,:) + tool.mvt_acc_MCS(IdxStep+1,:);
    tool.vel_MCS(IdxStep+1,:) = tool.vib_vel_MCS(IdxStep+1,:) + tool.mvt_vel_MCS(IdxStep+1,:);
    tool.pos_MCS(IdxStep+1,:) = tool.vib_pos_MCS(IdxStep+1,:) + tool.mvt_pos_MCS(IdxStep+1,:);
    
    % Projections to the Tool Coordinate System   
    tool.acc_TCS(IdxStep+1,:) = machine.ROT_MCS_TCS * tool.acc_MCS(IdxStep+1,:)';
    tool.vel_TCS(IdxStep+1,:) = machine.ROT_MCS_TCS * tool.vel_MCS(IdxStep+1,:)';
    tool.pos_TCS(IdxStep+1,:) = machine.ROT_MCS_TCS * tool.pos_MCS(IdxStep+1,:)';

    %% solve the vibration of the workpiece
    workpiece.vib_pos_MCS(IdxStep+1,1) = workpiece.vib_pos_MCS(IdxStep,1);
    workpiece.pos_MCS(IdxStep+1,2:3) = workpiece.vib_pos_MCS(IdxStep,2:3);
    
end % END of the main loop

figure(2)
subplot(2,1,2)
for IdxTooth=1:tool.Z % loop over all the teeth
    plot(solver.time,h_total(:,IdxTooth)*1000,'DisplayName',[ 'Tooth ', num2str(IdxTooth)])
    hold on;
end
xlabel('time (s)'); ylabel('chip thickness (mm)');
xlim([solver.time(1), solver.time(end)]); grid on;
