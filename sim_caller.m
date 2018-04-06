function sim_caller%(N85,N87,aprop,initModeVars)
% Interface for ARK45_GPE_sim

% File names
file_prefix     = 'ARK45_test';
groundstate_file= [file_prefix,'_groundstate'];
propagation_file= [file_prefix,'_propagation'];
pars_file       = [file_prefix,'_pars'];

% Set discretisation of spatial domain
n_x         = 2^6;
n_y         = 2^6;
n_z         = 2^6;

%% Experimental parameters
% parameters suffixed with _exp are the dimensional equivalents of
% non-dimensional (non-suffixed) parameters.

% Physical constants
hbar        = 1.0545718e-34; 
mass        = 1.4100e-25;       % Rb85, Rb87=1.4432e-25
gDQS        = 9.7959938810;     % At DQS, as measured by K. Hardman
bohr_radius = 5.2917721067e-11;

% Atom number
N           = 3.5e4;

% Set trap parameters
wx_exp  = 70*2*pi;
wy_exp  = 70*2*pi;
wz_exp  = 7*2*pi;

% Set size of spatial domain
size_x_exp  = 6e-6;
size_y_exp  = 6e-6;
size_z_exp  = 50e-6;

% Set interaction strengths
a_s         = 50*bohr_radius;            % scattering length
U_exp       = (4*pi*hbar^2*a_s/mass);
K3_im_exp   = -4.41e-41*hbar;   % 3-body recombination losses

% Set time stuff
t_max_exp           = 5e-3;            % total simulation time
dt_exp              = 0.05e-3;          % this is fixed if prop_mode == 'init'
sample_times_exp    = linspace(0,50e-3,100);% times to sample system

%% Nondimensionalise
Time    = 1/wx_exp;
Energy  = hbar/Time;    
Length  = sqrt(hbar*Time/mass);
KE_const= hbar^2 / (2*mass) / (Length^2 * Energy);
D_const = hbar / (4*mass) / (Length^2) * Time;

% Self-interaction
U           = U_exp/Length^3/Energy;

% 3 body losses
K3_im       = K3_im_exp/(Energy*Length^6);

% Trap stuff
wx              = wx_exp*Time;
wy              = wy_exp*Time;
wz              = wz_exp*Time;

% Set chemical potential
mu              = (15*N*U_exp*wx_exp*wy_exp*wz_exp/8/pi)^(2/5)*(mass/2)^(3/5)/Energy;

% Space stuff
size_x          = size_x_exp/Length;
size_y          = size_y_exp/Length;
size_z          = size_z_exp/Length;

% Time stuff
t_max           = t_max_exp/Time;        
dt              = dt_exp/Time;          
sample_times    = sample_times_exp/Time;

%% Define functions
% Trap
trap_fun    = @(x,y,z) 0.5*((wx*x).^2 + (wy*y).^2 + (wz*z).^2); % can also just supply array

% Scattering length (e.g. if it is changing over time)
f           = @(x) max(15*cos(1/(5e-3)/pi*x),10)*U*bohr_radius/a_s;

% Model function
model_fun   = @gpe_cq_3body_losses;
model_pars  = struct('K3_im',K3_im);

%% SET PARAMETERS
% Set trap geoemtry, number, and grid discretisation. trap_shift is shift
% along z-axis.
pars    = struct(...
            'N',N,...
            'trap_fun',trap_fun,...
            'size_x',size_x,...
            'size_y',size_y,...
            'size_z',size_z,...
            'n_x',n_x ,...
            'n_y',n_y ,...
            'n_z',n_z);
        
% Set normalisation parameters
pars    = appendfields(pars,...
            'KE_const',KE_const,...
            'D_const',D_const);
        
% Set model parameters
pars    = appendfields(pars,...
            'model_fun',model_fun,...
            'model_pars',model_pars);
        
% Set chemical potential (only used for initialisation)
pars    = appendfields(pars,...
            'chemical_potential',mu);
            
% Set ARK45 options
pars    = appendfields(pars,...
            'ark_tol',1e-1,...
            'ark_min_steps_before_change_dt',10,...
            'ark_min_step_size',1e-10,...
            'ark_method','DormandPrince');

% Set scattering lengths 
% note that U can be supplied as a constant or a function of the
% dimensionalised time.
pars    = appendfields(pars,...
            'U',f ...
            );

% Set three-body stuff
%3.3e-37
pars    = appendfields(pars,...
            'K3_im',K3_im);

% Set sim time and step size
pars    = appendfields(pars,...
            't_max',t_max,...
            'dt',dt,...
            'data_step_init',10,...
            'sample_times',sample_times);

% Set misc options
pars    = appendfields(pars,...
            'boundary_on',false,...
            'CUDA_on', false,...           
            'save_images_3D',false,...
            'save_images_2D',false,...
            'figures_on',true,...
            'write_video',false,...
            'debug_level',0);
        
        
% Set init/prop
pars    = appendfields(pars,...
            'prop_mode','init',...
            'init_file',groundstate_file,...
            'prop_file',propagation_file,...
            'overwriteInitFile',true,...
            'IP',true);
        
% Save parameters to file
save(pars_file,'pars')


% Initialise
% pars.prop_mode  = 'init';
% pars.U          = U;
% pars.t_max      = 20e-3/Time;
% ARK45_IP_GPE_3D_sim(pars)

% Propagate
pars.prop_mode  = 'prop';
pars.t_max      = 100e-3/Time;
pars.U          = @(t) U*max(cos(1/(5e-3)/pi*t*Time),0.5);
ARK45_GPE_3D_sim(pars)


end
