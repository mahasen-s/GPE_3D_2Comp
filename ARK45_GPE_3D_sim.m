function ARK45_GPE_3D_sim(pars)
% Should change the way arrays are written; should store piece wise on disk
% Switch for Cubic-Quintic GPE
% Functional form of integrator

% Cash-Carp and Fehlberg options
% Include checks for Levy-Friedrichs-Courant condition? ~dt/dx^2?

% Pass model as a reference
% Pass additional parameters as a struct

%% UNLOAD FROM PARAMETER STRUCTURE
% Font sizes for graphs
font_SmlSize= 10;
font_MedSize= 12;
font_LrgSize= 14;

% Normalisation parameters
KE_const        = pars.KE_const;
D_const         = pars.D_const;

% System parameters
N               = pars.N;
trap_fun        = pars.trap_fun;

% Model parameters
% model_fun       = pars.model_fun;
% model_pars      = pars.model_pars;

% Chemical potential
chemical_potential  = pars.chemical_potential;

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
ark_min_steps_before_change_dt  = pars.ark_min_steps_before_change_dt;
ark_min_step_size               = pars.ark_min_step_size;
ark_lo_ratio                    = 0.5;  % factor to reduce dt if step size too large
ark_hi_ratio                    = 1.9;  % factor to increase dt by if step too small, it's better if ark_lo_ratio*ark_hi_ratio != 1
ark_method                      = pars.ark_method;

% File IO
init_file       = pars.init_file;
prop_file       = pars.prop_file;

% Misc
prop_mode       = pars.prop_mode;
figures_on      = pars.figures_on;
boundary_on     = pars.boundary_on;
CUDA_on         = pars.CUDA_on;
save_images_3D  = pars.save_images_3D;
save_images_2D  = pars.save_images_2D;
write_video     = pars.write_video;
debug_level     = pars.debug_level;
%% Process inputs
% Set save file
switch prop_mode
    case 'init'
        save_file = init_file;
    case 'prop'
        save_file = prop_file;
end
% Check if plotting envrionment is available
if figures_on == true
    try
        main_fig = figure(1);
        clf
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
                if numel(U) == 1
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
disp_op     = get_disp_op;

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
        psi_fun     = real(sqrt((chemical_potential-potential)/abs(U)));
        
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
t_now    = 0;
if save_images_3D == true
    n_image              = numel(sampleTimes) + 1;
    image_array          = zeros(n_x,n_y,n_z,n_image) ;
    image_array(:,:,:,1) = psi_fun;
end

if save_images_2D == true
    n_image                 = numel(sampleTimes) + 1;
    image_array_xy          = zeros(n_x,n_y,n_image);
    image_array_yz          = zeros(n_y,n_z,n_image);
    image_array_zx          = zeros(n_z,n_x,n_image);
    
    [p_xy,p_yz,p_zx]        = get_proj(psi_fun);
    
    image_array_xy(:,:,1)   = p_xy;
    image_array_yz(:,:,1)   = p_yz;
    image_array_zx(:,:,1)   = p_zx;
end

% Atom number
N_now                   = get_N;
N_list                  = N_now;
mu_list                 = get_mu;

% U
if isequal(prop_mode,'prop')
    U_list  = U(0);
end

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

subtightplot(2,4,1,stp);
plot_proj_xy            = imagesc(x,y,p_xy);
xlabel('\boldmath$x$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$y$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
set(gca,'ydir','normal')
title('\boldmath$xy$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');

subtightplot(2,4,2,stp);
plot_proj_yz            = imagesc(y,z,p_yz);
xlabel('\boldmath$y$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$z$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
set(gca,'ydir','normal')
title('\boldmath$yz$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');

subtightplot(2,4,3,stp);
plot_proj_zx            = imagesc(x,z,p_zx);
xlabel('\boldmath$x$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$z$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
set(gca,'ydir','normal')
title('\boldmath$xz$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');

% Plot atom number
subtightplot(2,4,4,stp)
plot_N      = plot(time_vec,mu_list,'linewidth',1);%,'marker','o');
box on
grid on
title('\boldmath$N(t)$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
xlabel('\boldmath$t$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$N$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
xlim([0,t_max])

% Plot chemical potential
subtightplot(2,4,5,stp)
plot_mu     = plot(time_vec,mu_list,'linewidth',1);%,'marker','o');
box on
grid on
title('\boldmath$\mu(t)$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
xlabel('\boldmath$t$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
ylabel('\boldmath$\mu$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
xlim([0,t_max])

if isequal(prop_mode,'prop')
    % Plot U
    subtightplot(2,4,6,stp)
    plot_U     = plot(time_vec,U_list,'linewidth',1);
    box on
    grid on
    title('\boldmath$U(t)$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\boldmath$t$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\boldmath$U$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlim([0,t_max])
    
    % Plot dt
    subtightplot(2,4,7,stp)
    plot_dt    = plot(dt_time_list,dt_list,'linewidth',1);%,'marker','o');
    box on
    grid on
    title('\boldmath$dt(t)$','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\boldmath$t$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\boldmath$dt$','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlim([0,t_max])
    set(gca,'yscale','log')
end

%% PUSH ARRAYS TO GPU
if CUDA_on == true
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

% Prepare to capture video
if write_video == true
    set(main_fig,'color','w')
    vidObj  = VideoWriter([save_file,'_video.mj2'],'Archival');
end

% Pre-propagation
psi_fun         = fftn(psi_fun);

% Propagtion loop
tic
switch prop_mode
    case 'prop'
        ark_steps_without_change   = 0;
        t_next_samp_ind   = 0;    % index for next time to sample
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
                    updateSamplesAndPlot
                end
                
                % Advance to sample time
                propState
                get_sample_flag = false;
                
                % Get samples
                updateSamplesAndPlot
                
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
        U_now   = U;
        plotIter= 0;
        while t_now <= t_max  + dt
            
            % Do 1 step
            propStateRK4
            
            % Update t_now
            t_now        = t_now + abs(dt);
            
            % Update plot
            plotIter    = plotIter + 1;
            if mod(plotIter,data_step_init) == 0
                updateSamplesAndPlot
            end
            
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

% Close Video
if write_video == true
    close(vidObj);
end

% Post-propagation
psi_fun         = ifftn(psi_fun);

% Renormalise
if strcmp(prop_mode,'init') == true
    N_now  = get_N;
    psi_fun= psi_fun* sqrt(N/N_now);
end

% Gather
if CUDA_on == true
    psi_fun= gather(psi_fun);
end

% it's done
time_taken = toc;
fprintf('\nDone. Time taken \t=%4.3f\n',time_taken)

%% SAVE
fprintf('Saving  . . .\n')
switch prop_mode
    case 'init'
        save(save_file,...
            'time_vec',...
            'mu_list',...
            'pars',...
            'psi_fun');
    case 'prop'
        save(save_file,...
            'time_vec',...
            'N_list',...
            'U_list',...
            'dt_list',...
            'dt_time_list',...
            'pars',...
            'psi_fun');
end

if save_images_2D == true
    save(save_file,...
        'image_array_xy',...
        'image_array_yz',...
        'image_array_zx',...
        '-append')
end

if save_images_3D == true
    save(save_file,...
        'image_array',...
        '-append')
end

fprintf('Done\n\n')

%% GPE FUNCTION
    function out = gpeFun(psi_fun)
        dens    = abs(psi_fun).^2;
        out     = 1i*((-U_now*dens-K3_im*dens.^2-potential).*psi_fun);
    end

%% SUBFUNCTIONS
    function N  = get_N
        N   = trapz(x,trapz(y,trapz(z,abs(psi_fun).^2,3),2),1);
        if CUDA_on==true
            N   = gather(N);
        end
    end

    function mu  = get_mu
       psi_FFT = ( dx*dy*dz/((2*pi)^(3/2)) )*fftshift(fftn(fftshift(psi_fun)));
       KE = KE_const *trapz(kvec_x,trapz(kvec_y,trapz(kvec_z,(ifftshift(Kx).^2+ifftshift(Ky).^2+ifftshift(Kz).^2) .* abs(psi_FFT).^2,3),2),1);
       switch prop_mode
           case 'init'
               PE   = trapz(x,trapz(y,trapz(z,potential .* abs(psi_fun).^2 + U*abs(psi_fun).^4,3),2),1);
           case 'prop'
               PE   = trapz(x,trapz(y,trapz(z,potential .* abs(psi_fun).^2 + U(t_now)*abs(psi_fun).^4,3),2),1);
       end
       Ntot = trapz(x,trapz(y,trapz(z,abs(psi_fun).^2,3),2),1);
       mu = (KE + PE) / Ntot;
       if CUDA_on==true
           mu   = gather(mu);
       end
    end

    function disp_op = get_disp_op
        disp_op      = exp(-1i*D_const*dt*(Kx.^2+Ky.^2+Kz.^2));
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
            debugInfo('dt is approaching machine precision, integrator results may be unreliable',2);
        end
        
        % Add point at end of jump
        dt_list(end+1)       = new_dt;
        dt_time_list(end+1)  = t_now;
        
        % Recompute dispersion operator
        disp_op         = get_disp_op;
    end

    function [widths,coms] = getMoments
        % Get com
        com_x    = trapz(x,(x).*trapz(z,trapz(y,abs(ps2_85_data).^2,1),3),2)/N_now;
        com_y    = trapz(y,(y)'.*trapz(x,trapz(z,abs(ps2_85_data).^2,3),2),1)/N_now;
        z_temp      = reshape(z,[1,1,numel(z)]);
        com_z    = trapz(z,z_temp.*trapz(y,trapz(x,abs(ps2_85_data).^2,2),1),3)/N_now;
        
        % Get widths
        width_x  = (1/N_now)*trapz(x,(x.^2).*trapz(y,trapz(z,abs(ps2_85_data).^2,3),2),1);
        width_y  = (1/N_now)*trapz(y,(y.^2)'.*trapz(z,trapz(x,abs(ps2_85_data).^2,1),3),2);
        z_temp      = reshape(z,[1,1,numel(z)]);
        width_z  = (1/N_now)*trapz(z,(z_temp.^2).*trapz(x,trapz(y,abs(ps2_85_data).^2,2),1),3);
        
        % Store
        widths      = [ width_x;
            width_y;
            width_z];
        coms        = [ com_x;
            com_y;
            com_z];
    end

%% UPDATE PLOTS
    function updateSamplesAndPlot
        
        % UnFFT
        psi_fun = ifftn(psi_fun);
        
        % Update projections
        [p_xy,p_yz,p_zx]   = get_proj(psi_fun);
        
        if save_images_3D == true || save_images_2D==true
            sliceCounter     = sliceCounter+1;
        end
        
        if save_images_3D == true
            image_array(:,:,:,sliceCounter)  = psi_fun;
        end
        
        if save_images_2D == true
            image_array_xy(:,:,sliceCounter)    = p_xy;
            image_array_yz(:,:,sliceCounter)    = p_yz;
            image_array_zx(:,:,sliceCounter)    = p_zx;
        end
        
        % UPDATE QUANTITIES
        % Update time
        time_vec(end+1) = t_now;
        
        % Update particle counts
        N_now           = get_N;
        N_list(end+1)   = N_now;
        mu_list(end+1)  = get_mu;
        
        % Update U
        if isequal(prop_mode,'prop')
            U_list(end+1)   = U(t_now);
        end
        
        
        % UPDATE PLOTS
        if figures_on == true
            % Update prop plots
            if isequal(prop_mode,'prop') == true
                % Update dt
                set(plot_dt,'XData',dt_time_list,'YData',dt_list);
                
                % Update U
                set(plot_U,'XData',time_vec,'YData',U_list);
            end
            
            % Update mu, N
            set(plot_N,'XData',time_vec,'YData',N_list);
            set(plot_mu,'XData',time_vec,'YData',mu_list);
            
            % Update projections
            set(plot_proj_xy,'CData',p_xy);
            set(plot_proj_yz,'CData',p_yz);
            set(plot_proj_zx,'CData',p_zx);
        end
        
        % FFT
        psi_fun     = fftn(psi_fun);
        
        drawnow
        pause(0.001);
        
        % Write video frame
        if write_video == true
            frame   = getframe(main_fig);
            writeVideo(vidObj,frame)
        end
    end

%% PROP FUN
    function propState
        % Get current U
        U_now   = U(t_now);
        
        % Loop while step is not successful
        try_again   = true;     % try_again=false when step is successful
        
        % Store old state in case current step fails
        psi_fun_old = psi_fun;
        
        while try_again==true
            % Restore old values
            psi_fun     = psi_fun_old;
            
            % 1st half of dispersion
            psi_fun     = disp_op.*psi_fun;
            psi_fun     = ifftn(psi_fun);
            
            % Compute interaction using ARK45(DP)
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
            P_err = sum(abs(P_G(:)-P_H(:)))/numel(P_G);
            s     = min((dt*ark_tol./(2*(t_max)*P_err)).^(1/4));

            if s<1 || isnan(s)==true || isinf(s)==true
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
                    else
                        debugInfo(sprintf('Increasing step, dt->%4.3e',dtNew),2)
                    end
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
                    debugInfo(sprintf('Keeping step, t->%4.3e',t_now),4)
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
        psi_fun     = fftn(psi_fun);
        psi_fun     = disp_op.*psi_fun;
        
        
        % Update dt
        if dt_changed == true
            updateDt(dtNew)
        end
    end

    function propStateRK4
        % Apply 1st half of dispersion
        psi_fun     = disp_op.*psi_fun;
        psi_fun     = ifftn(psi_fun);
        
        % Interaction (RK4)
        P_A     = dt*gpeFun(psi_fun);
        P_B     = dt*gpeFun(psi_fun+ 1/2*P_A);
        P_C     = dt*gpeFun(psi_fun+ 1/2*P_B);
        P_D     = dt*gpeFun(psi_fun+ P_C);
        
        psi_fun = psi_fun + 1/6*(P_A + 2*P_B + 2*P_C + P_D);
        
        % Renormalise
        psi_fun = psi_fun*sqrt(N/get_N);
        
        % Apply 2nd half of dispersion
        psi_fun     = fftn(psi_fun);
        psi_fun     = disp_op.*psi_fun;
    end

%% MISC FUN
    function debugInfo(str,level)
        if debug_level>level
            fprintf([str,'\n'])
        end
    end
end

%% EXTERNAL FUNCTIONS
function mask = boundaries(n_i,n_j,n_k)
% probably should use something less bad
x = linspace(-1.1,1.1,n_i);
y = linspace(-1.1,1.1,n_j);
z = linspace(-1.1,1.1,n_k);
[X,Y,Z] = meshgrid(x,y,z);
mask = exp(-(X.^40+Y.^40+Z.^40));
end
