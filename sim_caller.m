function sim_caller%(N85,N87,aprop,initModeVars)
% Interface for ARK45_GPE_sim

% File names
file_prefix     = 'ARK45_test';
groundstate     = [file_prefix,'_groundstate'];
data_out        = [file_prefix,'_data'];
pars_out        = [file_prefix,'_pars'];
overwrite_flag  = true;

% Set system parameters
wx      = 70*2*pi;
wy      = 70*2*pi;
wz      = 7*2*pi;
scat_init = 50; % in units of a0

% Set constants
hbar        = 1.0545718e-34;
m85         = 1.4100e-25;
gDQS        = 9.7959938810;
bohr_radius = 5.2917721067e-11;

% Nondimensionalise
Time    = 1/omega_x(1);
Energy  = hbar/Time;    
Length  = sqrt(hbar*Time/m85);

% 85 self-interaction
a85         = scat_init*bohr_radius;
UaUnit85    = 4*pi*hbar^2*a85/m85;          % Units of interaction parameter
Ua85        = UaUnit85/Length^3/Energy*(scat_init/scat_init);

% 3 body losses
G3_im       = K3_im_prop*hbar/(Energy*Length^6);

% Define trap function
trap_fun    = @(x,y,z) 0.5*((wz*Time*x).^2 + (wy*Time*y).^2 + (wz*Time*z).^2);

%% SET PARAMETERS
% Set trap geoemtry, number, and grid discretisation. trap_shift is shift
% along z-axis.
pars    = struct(...
            'N',2500,...
            'trap_fun',trap_fun,...
            'size_x',6e-6,...
            'size_y',6e-6,...
            'size_z',50e-6,...
            'n_x',2^6 ,...
            'n_y',2^6 ,...
            'n_z',2^6,...
            'time_unit',1/wx);
            
% Set ARK45 options
pars    = appendfields(pars,...
            'tol',1e-6,...
            'ark_min_steps_before_change_dt',1);

% Set scattering lengths 
% note that scat_prop can be supplied as a constant or a function of the
% dimensionalised time.
pars    = appendfields(pars,...
            'scat_init',19,...
            'scat_prop',@(x) -50/(50e-3)^2*x.^2 + 50 ... 
            );

% Set three-body stuff
%3.3e-37
pars    = appendfields(pars,...
            'K3_re_init',0,...
            'K3_re_prop',0,...
            'K3_im_prop',0);%-4.41e-41);

% Set sim time and step size
pars    = appendfields(pars,...
            'Tmax_real_init',1e-3,...
            'dt_real_init',0.05e-5,...
            'data_step_init',100,...
            'Tmax_real_prop',10e-3,...
            'dt_real_prop',0.005e-5,...
            'sampleTimes',1e-3*linspace(0,50,1000),...
            'data_step_prop',1);

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
            'initMode','shift',...
            'initModeVars',[-5,5],...
            'init_file',groundstate,...
            'prop_file',data_out,...
            'overwriteInitFile',true);
        
% Physical constants
pars    = appendfields(pars,...
            'hbar',1.0545718e-34,...
            'm85',1.4100e-25,...
            'm87',1.4432e-25,...
            'gDQS',9.7959938810,...
            'bohr_radius',5.2917721067e-11);
        
% Set noise properties
pars    = appendfields(pars,...
            'noiseFrac',0.1,...
            'noiseStdDevs',[4e-7,4e-7,4e-7]);

        
% Save parameters to file
save(pars_out,'pars')

% a_85 s.t. 85-85 interaction is as strong as 85-87
a85_eqIntStr      = pars.m85/pars.m87*pars.scat_init_87;


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
