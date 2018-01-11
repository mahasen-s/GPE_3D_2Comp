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
N       = 2.5e4;

% Set trap parameters
wx_exp  = 70*2*pi;
wy_exp  = 70*2*pi;
wz_exp  = 7*2*pi;

% Set size of spatial domain
size_x_exp  = 6e-6;
size_y_exp  = 6e-6;
size_z_exp  = 50e-6;

% Set interaction strengths
a_s         = 50*a0;            % scattering length
U_exp       = (4*pi*hbar^2*a_s/mass);
K3_im_exp   = -4.41e-41*hbar;   % 3-body recombination losses

% Set time stuff
t_max_exp           = 50e-3;            % total simulation time
dt_exp              = 0.05e-5;          % this is fixed if prop_mode == 'init'
sample_times_exp    = linspace(0,50e-3);% times to sample system

%% Nondimensionalise
Time    = 1/wx_exp;
Energy  = hbar/Time;    
Length  = sqrt(hbar*Time/mass);

% Self-interaction
U           = U_exp/Length^3/Energy;

% 3 body losses
K3_im       = K3_im_exp/(Energy*Length^6);

% Trap stuff
wx              = wx_exp*Time;
wy              = wy_exp*Time;
wz              = wz_exp*Time;

% Space stuff
size_x_         = size_x_exp/Length;
size_y          = size_y_exp/Length;
size_z          = size_z_exp/Length;

% Time stuff
t_max           = Tmax_exp/Time;        
dt_exp          = dt_exp/Time;          
sample_times    = sample_times_exp/Time;

%% Define functions
% Trap
trap_fun    = @(x,y,z) 0.5*((wx*x).^2 + (wy*y).^2 + (wz*z).^2);

% Scattering length (e.g. if it is changing over time)
f           = @(x) max(15*cos(1/(5e-3)/pi*x),10)*U*a0/a_s;

%% SET PARAMETERS
% Set trap geoemtry, number, and grid discretisation. trap_shift is shift
% along z-axis.
pars    = struct(...
            'N',2500,...
            'trap_fun',trap_fun,...
            'size_x',size_x,...
            'size_y',size_y,...
            'size_z',size_z,...
            'n_x',n_x ,...
            'n_y',n_y ,...
            'n_z',n_z);
            
% Set ARK45 options
pars    = appendfields(pars,...
            'tol',1e-6,...
            'ark_min_steps_before_change_dt',1,...
            'ark_min_step_size',0,...
            'ark_method','DormandPrince');

% Set scattering lengths 
% note that scat_prop can be supplied as a constant or a function of the
% dimensionalised time.
pars    = appendfields(pars,...
            'a_s',f ...
            );

% Set three-body stuff
%3.3e-37
pars    = appendfields(pars,...
            'K3_im',K3_im);

% Set sim time and step size
pars    = appendfields(pars,...
            't_max',1e-3,...
            'dt',0.05e-5,...
            'data_step_init',100,...
            'sample_times',sample_times);

% Set misc options
pars    = appendfields(pars,...
            'bdd_on',false,...
            'CUDA_flag', true,...
            'useLogScale',false,...           
            'saveImages3D',false,...
            'saveImages2D',false,...
            'figuresOn',true,...
            'writeVideo',false);
        
        
% Set init/prop
pars    = appendfields(pars,...
            'propMode','init',...
            'init_file',groundstate_file,...
            'prop_file',propagation_file,...
            'overwriteInitFile',true);
        
% Save parameters to file
save(pars_file,'pars')

% Get groundstate
% gsFile              =  strrep([file_prefix,sprintf('8587_groundstate_as=%4.1fa0',a85_eqIntStr)]'.','p');
% pars.fileIn         = gsFile;
% CQSSFM_8587_3D_sim(pars)
% gsFile              =  strrep([file_prefix,sprintf('8587_groundstate_as=%4.1fa0_85Low',pars.scat_init_85)],'.','p');
% pars.fileIn         = gsFile;
% pars.initModeVars   = [1,-1]*10;
% CQSSFM_8587_3D_sim(pars)
% 
% gsFile              =  strrep([file_prefix,sprintf('8587_groundstate_as=%4.1fa0_87Low',pars.scat_init_85)],'.','p');
% pars.fileIn         = gsFile;
% pars.initModeVars   = [-1,1]*10;
% CQSSFM_8587_3D_sim(pars)

gsFile              =  strrep([file_prefix,sprintf('8587_groundstate_as=%4.1fa0_Sym',pars.scat_init_85)],'.','p');
pars.fileIn         = gsFile;
pars.initModeVars   = [0,0,1];
ARK45_SSFM_8587_3D_sim(pars)

% pars.propMode       = 'prop';
% ARK45_SSFM_8587_3D_sim(pars)

% gsFile              =  strrep([file_prefix,sprintf('8587_groundstate_as=%4.1fa0_Shift1',pars.scat_init_85)],'.','p');
% pars.fileIn         = gsFile;
% pars.initModeVars   = [-6,6,1];
% CQSSFM_8587_3D_sim(pars)
% 
% gsFile              =  strrep([file_prefix,sprintf('8587_groundstate_as=%4.1fa0_Shift2',pars.scat_init_85)],'.','p');
% pars.fileIn         = gsFile;
% pars.initModeVars   = [5,-5,1];
% CQSSFM_8587_3D_sim(pars)

% Prop
% parsIn                  =   pars;
% propFile                =  [file_prefix,sprintf('8587_waveguideProp_as=%4.1fa0',parsIn.scat_prop_85)];
% propFile                = strrep(propFile,'.','p');
% fprintf('Processing file %s . . . \n',propFile)
% parsIn.propMode         = 'prop';
% parsIn.fileOut          = propFile;
% CQSSFM_8587_3D_sim(parsIn)



end
