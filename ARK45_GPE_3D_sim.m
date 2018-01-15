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
ark_min_steps_before_change_dt  = pars.ark_ark_min_steps_before_change_dt;
ark_min_step_size               = pars.ark_ark_min_step_size;
ark_lo_ratio                    = 0.5;  % factor to reduce dt if step size too large
ark_hi_ratio                    = 1.9;  % factor to increase dt by if step too small, it's better if ark_lo_ratio*ark_hi_ratio != 1
ark_method                      = pars.ark_method;

% Misc
prop_mode       = pars.prop_mode;
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
        fprintf('Failed to create figures. Setting figures_on = false . . .\n')
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

% Check that U is correct
switch prop_mode
    case 'init'
        % U must be a numeric scalar
        if isnumeric(U)~=true || numel(U)~=1
            error('U must be a numeric scalar for prop_mode==init')
        end
    case 'prop'
        % If U is a numeric scalar, convert it to a function
        switch isnumeric(U)
            case true
                if numel(U) == 1 && pr
                    U = @(t) U;
                else
                    error('U must be a numeric scalar or a function a_s = f(t)')
                end
            otherwise
                if isa(U,'function_handle')==false
                    error('U must be a numeric scalar or a function a_s = f(t)')
                end
        end
end

% Switch to imaginary time for 'init', set K3_im to 0
switch prop_mode
    case 'init'
        dt      = -1i*dt;
        K3_im   = 0;
end

%% Define gpeFun
    function out = gpeFun(psi_fun)
        dens    = abs(psi_fun).^2;
        out     = 1i*((-U_now*dens-K3_im*dens.^2-potential).*psi_fun);
    end


%% Construct meshes
% Construct spatial vectors
dx      = 2*size_x/n_x;
dy      = 2*size_y/n_y;
dz      = 2*size_z/n_z;

x       = (-(n_x/2)*dx):dx:((n_x/2)*dx - dx);
y       = (-(n_y/2)*dy):dy:((n_y/2)*dy - dy);
z       = (-(n_z/2)*dz):dz:((n_z/2)*dz - dz);

[X,Y,Z] = ndgrid(x,y,z);

% Set up k-space coordinates
kvec_x      = pi/size_x*linspace(-n_x/2,n_x/2-1,n_x)';
kvec_y      = pi/size_y*linspace(-n_y/2,n_y/2-1,n_y)';
kvec_z      = pi/size_z*linspace(-n_z/2,n_z/2-1,n_z)';

% Create k-space meshes and centre zero-frequency component
[Kx,Ky,Kz]  = ndgrid(kvec_x,kvec_y,kvec_z);
Kx          = fftshift(Kx);
Ky          = fftshift(Ky);
Kz          = fftshift(Kz);

% Approximate dispersion in frequency space
disp_op     = exp(-1i*0.25*dt*(Kx.^2+Ky.^2+Kz.^2));

% Construct potential
    % TODO: Better abstraction for function input
switch class(trap_fun)
    case 'function_handle'
        try 
            potential   = trap_fun(X,Y,Z);
        catch me
            error('trap_fun improperly constructed')
        end
    case 'double'
        if isequal(size(X),size(trap_fun))
            potential = trap_fun;
            clear trap_fun
        else
            error('trap_fun is incorrectly sized')
        end
    otherwise
        error('trap_fun must be a function V = trap_fun(X,Y,Z) or a double array')
end

%% INITIALISE STATE
% Construct initial state, or load previous state
    % TODO: Better abstraction for function input
switch prop_mode
    case 'init'
    % Construct TF approximate ground state
    psi_fun     = real(sqrt((mu-potential)/abs(U(0))));
    
    case 'prop'
        % Load from file
        data    = load(init_file);
        psi_fun = data.psi_fun;
        
        clear data
end
        
% Clear coordinate meshes
clear X Y Z

%% SETUP ARRAYS FOR STORED QUANTITIES
% Create storage for images
if saveImages3D == true
   nImages              = numel(sampleTimes) + 1;
   imageArray           = zeros(n_x,n_y,n_z,nImages) ;
   imageArray(:,:,:,1)  = psi_fun;
end

if saveImages2D == true
    nImages                 = numel(sampleTimes) + 1;
    imageArray_xy           = zeros(n_x,n_y,nImages);
    imageArray_yz           = zeros(n_y,n_z,nImages);
    imageArray_zx           = zeros(n_z,n_x,nImages);   
    
    [p_xy,p_yz,p_zx]        = get_proj(psi_fun);
    
    imageArray_xy(:,:,1)    = p_xy;
    imageArray_yz(:,:,1)    = p_yz;
    imageArray_zx(:,:,1)    = p_zx;
end

% Atom number
N_now                   = 0;
N_list                  = getN;

% Time
time_vec                = 0;
dt_list                 = dt;
dt_time_list            = 0;

%% SETUP FIGURES
% subtightplot options
stp         = [0.07,0.08];
marg_h      = 1*[0.05,0.05];
marg_w      = [0.14,0.0];

% Plot projections
[p_xy,p_yz,p_zx]        = get_proj(psi_fun);

subtightplot(2,3,1);
plot_proj_xy            = imagesc(abs(p_xy));

subtightplot(2,3,2);
plot_proj_yz            = imagesc(abs(p_yz));

subtightplot(2,3,3);
plot_proj_zx            = imagesc(abs(p_zx));

% Plot atom number
subtightplot(2,3,4)
plot_N     = plot(time_vec,N_list,'linewidth',1);%,'marker','o');
box on 
grid on
title('\boldmath$dt(t)$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
xlabel('\boldmath$t$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$dt$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
xlim([0,t_max])

% Plot dt
plot_dt    = plot(dt_time_list,dt_list,'linewidth',1);%,'marker','o');
box on 
grid on
title('\boldmath$dt(t)$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
xlabel('\boldmath$t$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$dt$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
xlim([0,t_max])

%% PUSH ARRAYS TO GPU
if CUDA_flag == true
    x           = gpuArray(x);
    y           = gpuArray(y);
    z           = gpuArray(z);
    psi_fun     = gpuArray(psi_fun);
    disp_op     = gpuArray(disp_op);
    potential   = gpuArray(potential);
    if boundary_on  == true
        losses      = gpuArray(losses);    
    end
end

%% MAIN LOOP
fprintf('Progress: ');
prog_str = sprintf('00.00');
fprintf([prog_str,'%%']);

% Pre-propagation
psi_fun         = fftn(psi_fun);

% Propagtion loop
switch prop_mode
    case 'prop'
        t_next_samp_ind   = 0;    % index for next time to sample
        t_now    = 0;
        while t_now < t_max
            while t_next_samp_ind<numel(sample_times)
                % Get next sample time, tNext
                t_next_samp_ind = t_next_samp_ind + 1;
                t_next_samp     = sample_times(t_next_samp_ind);
                get_sample_flag = false; % flag to tell prop_state to sample on current step
                dt_is_shrunk_flag  = false; % flag to tell prop_state that dt was shrunk recently for sampling and that it should set dt = dt_pre_sample
                dt_pre_sample   = dt;
                
                % Propagate until sample time. At the end of this while loop,
                % tNow + dt = tNext
                while get_sample_flag == false
                    % Propagate until next step will put t_now at the correct
                    % sample time. propState will set get_sample_flag = true when
                    % this happens
                    propState
                end
                
                % Advance to sample time
                propState
                get_sample_flag = false;
                
                % Get samples
                updatePlots
                
                % Update progress
                prog_str = sprintf('%2.2f',t_now/t_max*100);
                if debug_level<1
                    % Update progress
                    fprintf(repmat('\b',1,length(prog_str)+1));
                    fprintf([prog_str,'%%']);
                else
                    fprintf(['Progress:\t',prog_str,'%%\n'])
                end
            end
        end
    case 'init'
        while t_now <= t_max  + dt
            
            % Do 1 step
            propState_RK4
            
            % Update tNow
            tNow        = tNow + abs(dt);
            
            % Update plot
            if mod(plotIter,data_step) == 0
                updatePlots
                
                % Write to video
                if write_video_flag == true && pars.figuresOn == true
                    frame = getframe(main_fig);
                    writeVideo(vidObj,frame);
                end
            end
            
            % Update progress
            fprintf(repmat('\b',1,length(str1)+1));
            str1 = sprintf('%2.2f',tNow/Tmax*100);
            fprintf([str1,'%%']);
            plotIter = plotIter+1;
            
            
        end
end



% Post-propagation
psi_fun         = ifftn(psi_fun);


%% PULL ARRAYS FROM GPU
psi_fun         = gather(psi_fun);


%% SAVE
save



%% Subfunctions
    function N  = getN
        N   = trapz(x,trapz(y,trapz(z,abs(psi_fun).^2,3),2),1);
        if CUDA_on==true
            N   = gather(N);
        end
        N_now   = N;
    end

    function [p_xy,p_yz,p_zx]   = get_proj(psi_fun)
        % Get 2D projections
        dens    = abs(psi_fun).^2;
        p_xy    = squeeze(sum(dens,3));
        p_yz    = squeeze(sum(dens,1));
        p_zx    = squeeze(sum(dens,2));
        
        if CUDA_on==true
            p_xy   = gather(p_xy);
            p_yz   = gather(p_yz);
            p_zx   = gather(p_zx);
        end
    end

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
                
            elseif s>2 && ark_steps_without_change >= ark_min_steps_before_change_dt
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
                    dtNew           = max(t_next_samp-t_now,ark_min_step_size);
                    
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
                    dtNew           = max(t_next_samp-t_now,ark_min_step_size);
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
    
    function updatePlots
        % Update things
        time_vec(end+1) = t_now;
        N_list(end+1)   = getN;
        
        % Update projections
        [p_xy,p_yz,p_zx]        = get_proj(psi_fun);
        set(plot_proj_xy,'CData',p_xy);
        set(plot_proj_yz,'CData',p_yz);
        set(plot_proj_zx,'CData',p_zx);
        
        % Update dt_list plot
        set(plot_dt,'XData',dt_time_list,'YData',dt_list)
        
        % Save images
        if saveImages3D==true || saveImages2D==true
            sliceCounter     = sliceCounter + 1;
        end
        if saveImages3D == true
            imageArray(:,:,:,t_next_samp_ind+1)  = psi_fun;
        end
        
        if saveImages2D == true
            imageArray_xy(:,:,t_next_samp_ind+1)    = p_xy;
            imageArray_yz(:,:,t_next_samp_ind+1)    = p_yz;
            imageArray_zx(:,:,t_next_samp_ind+1)    = p_zx;
        end
        
        % Pause biefly and draw
        drawnow
        pause(0.001)
    end

end
