function sim_caller_250517_fun%(N85,N87,aprop,initModeVars)
% MI in 85/87 mixture
% Start in miscible 85-87 regime. Switch to waveguide, jump scattering
% length to negative value. Does MI still occur? Does the 85 form a 1D
% lattice which traps 87?


% File names
file_prefix     = 'ARK45_test';
groundstate     = [file_prefix,'_groundstate'];
data_out        = [file_prefix,'_data'];
pars_out        = [file_prefix,'_pars'];
overwrite_flag  = true;

% Set constants


%% SET PARAMETERS
% Set trap geoemtry, number, and grid discretisation. trap_shift is shift
% along z-axis.
pars    = struct(...
            'N_85',2500,...
            'N_87',2500,...
            'freq_x',[70,70],...
            'freq_y',[70,70],...
            'freq_z',[7,7],...
            'freq_z_prop',[0,0],...
            'trap_shift', 0,...
            'size_x',6e-6,...
            'size_y',6e-6,...
            'size_z',50e-6,...
            'n_x',2^6 ,...
            'n_y',2^6 ,...
            'n_z',2^6 ,...
            'wgTilt',0,...%0.46*pi/180,...
            'gravityOn', true,...
            'tol',1e-6,...
            'ark_minStepsBeforeChangeDt',1);

% Set scattering lengths 
pars    = appendfields(pars,...
            'RampTime',0,...
            'scat_init_85',19,...
            'scat_init_87',108,...
            'scat_init_8587',214,...
            'scat_prop_85',19,...
            'scat_prop_87',108,...
            'scat_prop_8587',214,...
            'scat_exp_85',-200,...
            'scat_exp_87',100,...
            'scat_exp_8587',214);

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
            'sampleTimes',1e-3*linspace(0,10,1000),...
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
            'fileIn',groundstate,...
            'fileOut',data_out,...
            'overwriteFileIn',true);
        
% Set properties for pre-imaging expansion. expFactor_r, expFactor_z 
% indicate size of new array relative to original to be used for the
% expanded image.
pars    = appendfields(pars,...
            'ExpandOn',false,...
            'expFactor_r',3,...
            'expFactor_z',3,...
            'expTime',0e-3,...
            'expPostTime',0e-3,...
            'expPostStep',0e-3);
        
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
