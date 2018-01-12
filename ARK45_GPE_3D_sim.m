function ARK45_GPE_3D_sim(pars)
% Explicit visualisation controls
% Checks for whether figure was initialised done
% Video writing contingent on this check done

% Nondimensionalised interface done

% Switch for Cubic-Quintic GPE

% Functional form of integrator

% Cash-Carp and Fehlberg options

% trap can be provided either as a function or as an array

% Include checks for Levy-Friedrichs-Courant condition? ~dt/dx^2?

%% UNLOAD FROM PARAMETER STRUCTURE

% System parameters
N               = pars.N_85;
trap_fun        = pars.trap_fun;

% Interaction parameters
U               = pars.U;
K3_im           = pars.K3_im;

% Time
t_max           = pars.t_max;
dt              = pars.dt;
sample_times    = pars.sample_times;
data_step_init  = pars.data_step_init;

% Spatial size and discretisation
size_x          = pars.size_x;
size_y          = pars.size_y;
size_z          = pars.size_z;                     
n_x             = pars.n_x;
n_y             = pars.n_y;
n_z             = pars.n_z;

% ARK45 options
ark_tol                         = pars.ark_tol;
ark_min_steps_before_change_dt  = pars.ark_tol;
ark_min_step_size               = pars.ark_tol;
ark_method                      = pars.ark_tol;

% Misc
figures_on      = pars.figures_on;
boundary_on     = pars.boundary_on;
CUDA_on         = pars.CUDA_on;
use_log_scale   = pars.use_log_scale;
save_images_3D  = pars.save_images_3D;
save_images_2D  = pars.save_images_2D;
write_video     = pars.write_video;

%% Process inputs
% Check if plotting envrionment is available
if figures_on == true
    try
        h_fig = figure(1);
    catch me
        fprintf('Failed to create figures. Disabling visualisations . . .\n')
        figures_on = false;
    end
end

% If figures_on==false, disable all visualisation flags
if figures_on==false
    if save_images_2D==true || save_images_3D==true
        fprintf('Cannot save images since figures_on=false\n')
    end
    if write_video==true
        fprintf('Cannot write video since figures_on=false\n')
    end
end

% If a_s is a numeric scalar, convert it to a function
switch isnumeric(U)
    case true
        if numel(U) == 1
            U = @(t) U;
        else
            error('a_s must be a numeric scalar or a function a_s = f(t)')
        end
    otherwise
        if isa(U,'function_handle')==false
            error('a_s must be a numeric scalar or a function a_s = f(t)')
        end
end

%% Define gpeFun
    function out = gpeFun(psi_fun)
        dens    = abs(psi_fun).^2;
        out     = 1i*((-U_now*dens-K3_im*dens.^2-potential).*psi_fun);
    end


%% Subfunctions
	function updateDt(new_dt)
        % Update dt and dispersion operators
        
        % Add point at start of jump
        dt_list(end+1)       = dt;
        dt_time_list(end+1)  = t_now;
        
        % Update dt
        dt  = new_dt;
        
        % Check if smaller than min step
        if new_dt<16*eps(t_now)
            debugInfo('dt is approaching machine precision, integrator results may be unreliable',-1);
        end
        
        % Add point at end of jump
        dt_list(end+1)       = new_dt;
        dt_time_list(end+1)  = t_now;
        
        % Recompute dispersion operator
        disp_op         = exp(-1i*0.25*dt*(Kx.^2+Ky.^2+Kz.^2));
    end
    
    function propState
        % Loop while step is not successful
        try_again   = true;     % try_again=false when step is successful
        
        while try_again==true
            % Store old state in case current step fails
            psi_fun_old = psi_fun;
            
            % 1st half of dispersion
            psi_fun     = disp_op.*psi_fun;
            psi_fun     = ifftn(psi_fun);
            
            % Compute interaction using ARK45(DP)
            try_again   = true;     % try_again==false when step is successful
            dt_changed  = false;    % step successful, but dt needs to be increased for efficiency, or reduced for sampling
            
            % ARK45(DP)
            P_A   = dt*gpeFun(psi_fun);
            P_B   = dt*gpeFun(psi_fun +          1/5*P_A);
            P_C   = dt*gpeFun(psi_fun +         3/40*P_A +       9/40*P_B);
            P_D   = dt*gpeFun(psi_fun +        44/45*P_A -      56/15*P_B +       32/9*P_C);
            P_E   = dt*gpeFun(psi_fun +   19372/6561*P_A - 25360/2187*P_B + 64448/6561*P_C - 212/729*P_D);
            P_F   = dt*gpeFun(psi_fun +    9017/3168*P_A -     355/33*P_B + 46732/5247*P_C +	49/176*P_D -    5103/18656*P_E);
            
            % 5th, 4th order estimates
            P_G   = dt*gpeFun(psi_fun +       35/384*P_A +          0*P_B +   500/1113*P_C + 125/192*P_D -     2187/6784*P_E +    11/84*P_F           );
            P_H   = dt*gpeFun(psi_fun +   5179/57600*P_A +          0*P_B + 7571/16695*P_C + 393/640*P_D -  92097/339200*P_E + 187/2100*P_F + 1/40*P_G);
            
            % Calculate error
            P_err = sum(abs(P_G(:)-P_H(:)));%/numel(P_G);
            s     = min((dt*ark_tol./(2*(t_max)*P_err)).^(1/5));
            
            if s<0
                % Previous step was too large. Halve and try again.
                updateDt(ark_lo_ratio*dt)
                try_again  = true;
                debugInfo(sprintf('Reducing step, P_err=%4.3e, s=%4.3e, dt=%4.3e',P_err,s,dt),2)
                
                % Reset counter
                ark_steps_without_change   = 0;
                
            elseif s>2 && ark_steps_without_change >= ark_steps_before_change_dt
                % Previous step was too small.
                % Update ps2 and increase time step-size for next step,
                % unless sampling
                dt_changed      = true;
                t_now            = t_now + dt;
                
                if t_next_samp-t_now < ark_hi_ratio*dt && get_sample_flag == false
                    % If we increase the time step too much, we would miss
                    % the sample time, so increase dt appropriately
                    % However, don't make the step smaller than min_step
                    % Store previous step size, and restore afterwards
                    % Also make sure that we aren't already sampling
                    
                    % Store pre-sample dt
                    dt_pre_sample   = dt;
                    
                    % Update dt, being careful to avoid min_step size
                    dtNew           = max(t_next_samp-t_now,min_step_size);
                    
                    % Update flags
                    dt_is_shrunk_flag = true;
                    get_sample_flag   = true;
                    
                    debugInfo(sprintf('Increasing step, dt->%4.3e, but sampling',dtNew),2)
                else
                    % Increase dt
                    dtNew  = ark_hi_ratio*dt;
                    
                    % If dt has been shrunk for sampling in the previous
                    % step, restore it to its original size
                    if dt_is_shrunk_flag==true
                        % Restore dt
                        dtNew   = dt_pre_sample;
                        
                        % Update flags
                        dt_is_shrunk_flag = false;
                        
                        debugInfo(sprintf('Restoring dt to pre-sample value, dt->%4.3e',dtNew),2)
                    end
                    
                    debugInfo(sprintf('Increasing step, dt->%4.3e',dtNew),2)
                end
                try_again  = false;
                
                % Reset counter
                ark_steps_without_change   = 0;
            else
                % Step-size was about right OR dt has been changed
                % recently
                % Update ps2 and leave dt unchanged, unless getting sample,
                % or restoring step size
                t_now                = t_now + dt;
                if t_next_samp-t_now < dt && get_sample_flag == false
                    % If we leave dt unchanged, we would miss
                    % the sample time, so reduce dt appropriately
                    % However, don't make the step smaller than min_step
                    % Also make sure that we aren't already sampling
                    dt_changed      = true;
                    
                    % Update dt, being careful to avoid min_step size
                    dtNew           = max(t_next_samp-t_now,min_step_size);
                    get_sample_flag   = true;
                    debugInfo(sprintf('Keeping step, but reducing dt for sampling, dt->%4.3e t->%4.3e',dtNew,t_now),2)
                    
                elseif dt_is_shrunk_flag==true
                    % Update dt
                    dtNew   = dt_pre_sample;
                    
                    % Update flags
                    dt_is_shrunk_flag   = false;
                    dt_changed          = true;
                    
                    debugInfo(sprintf('Restoring dt to pre-sample value, dt->%4.3e t->%4.3e',dtNew,t_now),2)
                else
                    debugInfo(sprintf('Keeping step, t->%4.3e',dtNew),4)
                end
                try_again  = false;
                
                % Increment unchanged counter
                ark_steps_without_change   = ark_steps_without_change + 1;
            end
            
            if try_again==false
                % Step was successful
                psi_fun   = psi_fun + P_G;
            end
        end
        % 2nd half of dispersion
        psi_fun     = disp_op.*psi_fun;
        psi_fun     = ifftn(psi_fun);
            
        % Update dt
        if dt_changed == true
            updateDt(dtNew)
        end
    end


end