function ARK45_SSFM_8587_3D_sim(pars)
%% Simulates interaction between an 85 and 87 condensate
% Accepts 'init','prop' as modes
% Callable function for propagating a groundstate using cylindrically 
% symmetric assumptions. Takes arguments *file_state_in* (name of file_out
% produced by groundstate code), *file_data_out* (name of output file 
% containing radial and longitudinal cloud widths and particle count) and 
% *pars* (structure containing all variables and parameters)

% Includes extra steps before imaging: free-space (a_s=0) expansion  
% (with g_3) followed by a period of post-expansion with a_s =Go 
% pars.scat_exp

% Shifts CoM to observe CoM Oscillations

% Are there loss terms for Rb87?

% Parameter for number of samples to take should be fixed
% need to include taking data at fixed time steps
% Add half sine ramp

% Put back time vector for updatedt Done
% set up figure for dt vs time Done
% fix figure arraay to accommodate new fiigures Done
% fix aprameters of sim_caller 
% make sure samppling works properly

%% UNLOAD FROM PARAMETER STRUCTURE

% Properties of trap, number, scattering length
N_85            = pars.N_85;                          
N_87            = pars.N_87;                         
freq_x          = pars.freq_x;
freq_y          = pars.freq_y;
freq_z          = pars.freq_z;                     
freq_z_prop     = pars.freq_z_prop;           
trap_shift      = pars.trap_shift;
size_x          = pars.size_x;
size_y          = pars.size_y;
size_z          = pars.size_z;                     
n_x             = pars.n_x;
n_y             = pars.n_y;
n_z             = pars.n_z; 

RampTime        = pars.RampTime;                 
scat_init_85    = pars.scat_init_85;         
scat_init_87    = pars.scat_init_87;         
scat_init_8587  = pars.scat_init_8587;     
scat_prop_85    = pars.scat_prop_85;         
scat_prop_87    = pars.scat_prop_87;         
scat_prop_8587  = pars.scat_prop_8587;     
scat_exp_85     = pars.scat_exp_85;           
scat_exp_87     = pars.scat_exp_87;           
scat_exp_8587   = pars.scat_exp_8587;  

K3_re_init      = pars.K3_re_init;             
K3_re_prop      = pars.K3_re_prop;             
K3_im_prop      = pars.K3_im_prop;    

Tmax_real_init  = pars.Tmax_real_init;     
dt_real_init    = pars.dt_real_init;         
Tmax_real_prop  = pars.Tmax_real_prop;     
dt_real_prop    = pars.dt_real_prop;

gravityOn       = pars.gravityOn;
sampleTimes     = pars.sampleTimes;
bdd_on          = pars.bdd_on;                     
CUDA_flag       = pars.CUDA_flag;               
useLogScale     = pars.useLogScale;           
saveImages3D    = pars.saveImages3D;
saveImages2D    = pars.saveImages2D;

% ARK vars
arkStepsWithoutChange = 0;

% Physical constants
hbar            = pars.hbar;
m85             = pars.m85;
m87             = pars.m87;
gDQS            = pars.gDQS;
bohr_radius     = pars.bohr_radius;    

% Set noise properties
noiseFrac       = pars.noiseFrac;       % Normalised size of added noise
noiseStdDevs    = pars.noiseStdDevs;    % Std devs in r and z (respectively) for gaussian modulation of noise spatial frequencies

% Waveguide tilt
wgTilt          = pars.wgTilt;

% Clean up sampleTimes (sort and remove zeros)
sampleTimes(sampleTimes==0) = [];
sampleTimes     = sort(sampleTimes);

% Determine if we are finding groundstate or propagating
propMode        = pars.propMode;
fileIn          = pars.fileIn;
fileOut         = pars.fileOut;
overwriteFileIn = pars.overwriteFileIn;

% Determine whether we are ramping to negative a or not
if pars.RampTime<= 0 || scat_init_85==scat_prop_85 || strcmp(propMode,'init')==true
    RampOn          = false;
else
    RampOn          = true;
    rampTime_real   = pars.RampTime;
end

% Set values depending on mode
switch propMode
    case 'init'
        data_step   = pars.data_step_init;
        save_file   = fileIn;
        if exist([save_file,'.mat'],'file')~=0 && overwriteFileIn == false
            disp('Groundstate already exists')
            return
        end
    case 'prop'
        data_step   = pars.data_step_prop;
        save_file   = fileOut;
end

% Capture video
if pars.writeVideo == true && strcmp(propMode,'prop')==1
   writeVideo_flag   = true; 
else
   writeVideo_flag   = false;
end
        

%% SET UP
 % Font sizes for graphs
font_SmlSize= 10;              
font_MedSize= 12;
font_LrgSize= 14;

% Set up boundary, if required
if bdd_on == true
    losses      = boundaries(n_x,n_y,n_z);
end

% Trap frequencies
omega_x     = 2*pi*freq_x;
omega_y     = 2*pi*freq_y;
switch propMode
    case 'init'
        omega_z     = 2*pi*freq_z;
    case 'prop'
        omega_z     = 2*pi*freq_z_prop;
end

% Nondimensionalised units
Time    = 1/omega_x(1);
Energy  = hbar/Time;    %Energy  = hbar*omega_x;
Length  = sqrt(hbar*Time/m85);

% Nondimensionalise time
switch propMode
    case 'init'
        Tmax    = Tmax_real_init/Time;
        dt      = -1i*dt_real_init/Time;
        Ttotal  = Tmax_real_init;
    case 'prop'
        Tmax    = Tmax_real_prop/Time;
        dt      = dt_real_prop/Time;
        Ttotal  = Tmax_real_prop;
end

% Nondimensionalise sizes
size_x  = size_x/Length;
size_y  = size_y/Length;
size_z  = size_z/Length;

% 85 self-interaction
a85         = scat_init_85*bohr_radius; % Bohr radius * scattering length    
UaUnit85    = 4*pi*hbar^2*a85/m85;          % Units of interaction parameter

% 87 self-interaction
a87         = scat_init_87*bohr_radius;     
Ua87        = (4*pi*hbar^2*a87/m87)/(Energy*Length^3); 

% Cross-interaction
a8587       = scat_init_8587*bohr_radius;     
m_eff       = 2/(1/m85+1/m87);                % Effective mass
Ua8587      = (4*pi*hbar^2*a8587/m_eff)/(Energy*Length^3); % Missing factor of 2? Or Paltins reduced mass is wrong?
 
% Set interaction and losses
switch propMode
    case 'init'
        G3_im       = 0;
        Ua85        = UaUnit85/Length^3/Energy*(scat_init_85/scat_init_85);
    case 'prop'
        G3_im       = K3_im_prop*hbar/(Energy*Length^6);
        Ua85        = UaUnit85/Length^3/Energy*(scat_prop_85/scat_init_85);
end
G3_re       = K3_re_prop*hbar/(Energy*Length^6);
G3          = G3_re + 1i*G3_im;

% Construct spatial vectors
dx      = 2*size_x/(n_x-1);
dy      = 2*size_y/(n_y-1);
dz      = 2*size_z/(n_z-1);
x       = -size_x:dx:size_x;
y       = -size_y:dy:size_y;
z       = -size_z:dz:size_z;

% Construct time vector
t       = 0:abs(dt):Tmax;
n_t     = numel(t);

% Construct potential for propagation
[X,Y,Z]                                         = meshgrid(x,y,z);
[potential85,potential87,trap85,trap87]         = getPotential(omega_x,omega_y,omega_z);
% figure(2)
% maxis = [inf,-inf];
% for i=1:3
%     subplot(2,3,i)
%     imagesc(squeeze(sum(potential85,i)))
%     cax = caxis;
%     if cax(1) < maxis(1)
%         maxis(1) = cax(1);
%     elseif cax(2) > maxis(2)
%         maxis(2) = cax(2);
%     end
%     
%     subplot(2,3,3+i)
%     imagesc(squeeze(sum(potential87,i)))
%     cax = caxis;
%     if cax(1) < maxis(1)
%         maxis(1) = cax(1);
%     elseif cax(2) > maxis(2)
%         maxis(2) = cax(2);
%     end
% end
% for i=1:6
%     subplot(2,3,i)
%     colormap(jet(2048))
%     caxis(maxis)
% end
% maxis
% pause

% Set up state
switch propMode
    case 'init'
        % Construct Thomas-Fermi approximate ground state
        UaUnitTF_85 = 4*pi*hbar^2*a85/m85;
        UaTF_85     = UaUnitTF_85/(Energy*Length^3);
        mu          = (15*N_85*UaUnitTF_85*omega_x(1)*omega_y(1)*omega_z(1)/8/pi)^(2/5)*(m85/2)^(3/5)/Energy;
        ps2_85      = real(sqrt((mu-trap85)/abs(UaTF_85)));
        
        UaUnitTF_87 = 4*pi*hbar^2*a87/m87;
        UaTF_87     = UaUnitTF_87/(Energy*Length^3);
        mu87        = (15*N_87*UaUnitTF_87*omega_x(2)*omega_y(2)*omega_z(2)/8/pi)^(2/5)*(m87/2)^(3/5)/Energy; %Should m85 be m87? Doing this ruins normalisation....
        
        switch pars.initMode
            % TF
            case 'tf'
            ps2_87      = real(sqrt((mu87-trap87)/abs(UaTF_87)));

            % Torus
            case 'torus'
                torus_rad   = 2.5e-6/Length;
                torus_wall  = 1e-6/Length;
                trap_t      = 0.5*((omega_z*Time*Z).^2);
                ps2_87_z    = real(sqrt((mu-m87/m85*trap_t)/abs(UaTF_87)));             % thomas-fermi in z
                ps2_87_xy   = exp(-(sqrt(X.^2+Y.^2)-torus_rad).^2/(2*torus_wall)^2);    % gaussian walled annulus in x,y
                ps2_87      = ps2_87_xy.*ps2_87_z;
                [~,n87_t]   = getN;
                ps2_87      = ps2_87*sqrt(N_87/n87_t);
%                 figure(2)
%                 clf
%                 subtightplot(1,3,1)
%                 imagesc(abs(squeeze(sum(trap_t,2))));
%                 subtightplot(1,3,2)
%                 imagesc(abs(squeeze(sum(ps2_87_z))));
%                 subtightplot(1,3,3)
%                 imagesc(abs(squeeze(sum(ps2_87_xy))));
                
                clear trap_t ps2_87_z ps2_87_xy

            % Two gaussians
            case 'gaussian'
                ps2_87      =   (exp(-(Z-20e-6/Length).^2/(2*(10e-6/Length)^2)) + exp(-(Z+20e-6/Length).^2/(2*(10e-6/Length)^2))).*...
                    exp(-(X.^2+Y.^2)/(2*(2e-6/Length)^2));
                
            case 'shift'
                ps2_87 = circshift(real(sqrt((mu87-trap87)/abs(UaTF_87))),pars.initModeVars(1),pars.initModeVars(3));
                ps2_85 = circshift(ps2_85,pars.initModeVars(2),pars.initModeVars(3));
        end
    case 'prop'
        % Load initial state
        dataIn  = load(fileIn);
        ps2_85  = dataIn.ps2_85;
        ps2_87  = dataIn.ps2_87;
        
        % Clear
        clear dataIn
end

% potential85 = trap85;
% potential87 = trap87;
clear X Y Z trap85 trap87

%% Set up k-space coordinates
% Set up k-space coordinates
kvec_x      = pi/size_x*linspace(-n_x/2,n_x/2-1,n_x)';
kvec_y      = pi/size_y*linspace(-n_y/2,n_y/2-1,n_y)';
kvec_z      = pi/size_z*linspace(-n_z/2,n_z/2-1,n_z)';

% Create k-space meshes and centre zero-frequency component
[Kx,Ky,Kz]  = meshgrid(kvec_x,kvec_y,kvec_z);
Kx          = fftshift(Kx);
Ky          = fftshift(Ky);
Kz          = fftshift(Kz);


%% ADD NOISE
if strcmp(propMode,'prop') == 1
    % Create noise profile
    noise_sig_x     = noiseStdDevs(1)/dx/Length;
    noise_sig_y     = noiseStdDevs(2)/dy/Length;
    noise_sig_z     = noiseStdDevs(3)/dz/Length;
    noiseProfile    = 2*(1-rand(size(ps2_85))) + 2*1i*(1-rand(size(ps2_85)));
    noiseProfile    = fftshift(fftshift(fft2(noiseProfile)).*...
                        exp(- Kx.^2/(2*noise_sig_x.^2)...
                            - Ky.^2/(2*noise_sig_y.^2)...
                            - Kz.^2/(2*noise_sig_z.^2)));
    noiseProfile    = ifft2(noiseProfile);

    % Add to psi
    ps2_85_orig     = ps2_85;
    ps2_85          = ps2_85.*(1+noiseFrac*noiseProfile);
    
    % Clear vars
    clear noiseProfile

    % Renormalise
    [n_85_temp,~]   = getN;
    ps2_85          = ps2_85 * sqrt(N_85/n_85_temp); 
end

%% SETUP EXPANSION STUFF
if pars.ExpandOn == true
    
    % Expansion Ua unit
    Ua_exp_85  = UaUnit85/Length^3/Energy*(scat_exp_85/scat_init_85);
    
    % Get expansion time
    dt_PostExp  = dt_PostExp_real/Time;
    expTime = expTime_real/Time;
    expPostTime = expPostTime_real/Time;
    
    % Get expansion time vector
    t_exp   = 0:abs(dt):expTime;
    t_postexp = 0:dt_PostExp:expPostTime;
    n_t_exp = numel(t_exp);
    n_t_postexp = numel(t_postexp);
end

%% SETUP RAMP STUFF
if RampOn == true
    rampTime    = rampTime_real/Time;
    rampUa_fun  = @(x) UaUnit85/Length^3/Energy*((scat_init_85 + x*(scat_prop_85-scat_init_85)/rampTime)/scat_init_85);
end

%% PREPARE FOR PROPAGATION
% Approximate dispersion in frequency space
disp_op85       = exp(-1i*0.25*dt*(Kx.^2+Ky.^2+Kz.^2));
disp_op87       = exp(-1i*0.25*dt*(Kx.^2+Ky.^2+Kz.^2)*m85/m87);

if pars.ExpandOn == true
    disp_op_PostExp     = exp(-1i*(Kr.^2+Kz.^2)*dt_PostExp*0.25);
end

% Clear K
if strcmp(propMode,'prop') == false
    clear Kx Ky Kz
end

%% SETUP MONITORED QUANTITIES
main_fig = figure(1);
clf
set(main_fig,'renderer','opengl','color','w')

% Get data of unexpanded cloud
ps2_85_data                 = ps2_85;
ps2_87_data                 = ps2_87;

[n85,n87]                   = getN; % n85, n87 are current values of N; N_85, N_87 are initial
Nlist_85                    = n85;
Nlist_87                    = n87;

[E_85,E_87]                 = getE;
Elist_85                    = E_85;%NaN;
Elist_87                    = E_87;%NaN;

[widths,coms]               = getMoments;
width_x_list_85             = widths(1,1);
width_y_list_85             = widths(2,1);
width_z_list_85             = widths(3,1);
com_x_list_85               = coms(1,1);
com_y_list_85               = coms(2,1);
com_z_list_85               = coms(3,1);

width_x_list_87             = widths(1,2);
width_y_list_87             = widths(2,2);
width_z_list_87             = widths(3,2);
com_x_list_87               = coms(1,2);
com_y_list_87               = coms(2,2);
com_z_list_87               = coms(3,2);

time_vec                    = 0;
dtList                      = dt;
dtTimeList                  = 0;

% Get data of expanded cloud
if pars.ExpandOn == true
    % FFT and Hankel
    ps2_85_exp     = HankelFFT(ps2_85_exp);
    ps2_87_exp     = HankelFFT(ps2_87_exp);

    % Expand
    if pars.ExpandOn == true
        expandCloud
    end

    % UnHankel and UnFFT
    ps2_85_exp     = iHankelFFT(ps2_85_exp);
    ps2_87_exp     = iHankelFFT(ps2_87_exp);
    if CUDA_flag == true
        ps2_85_exp= gather(ps2_85_exp);
        ps2_87_exp= gather(ps2_87_exp);
    end

    [ne85,ne87]         = getExpandedN;
    Nlist_exp_85        = ne85;
    Nlist_exp_87        = ne87;
    [width_z_exp, width_r_exp,com_z_exp]  = getExpandedWidths85;
    width_z_exp_list_85         = width_z_exp;
    width_r_exp_list_85         = width_r_exp;
    com_z_exp_list_85           = com_z_exp;
    [width_z_exp, width_r_exp,com_z_exp]  = getExpandedWidths87;
    width_z_exp_list_87         = width_z_exp;
    width_r_exp_list_87         = width_r_exp;
    com_z_exp_list_87           = com_z_exp;
end

% Create storage for images
if saveImages3D == true && strcmp(propMode,'init') == false
   nImages            = numel(sampleTimes) + 1;
   imageArray85       = zeros(n_x,n_y,n_z,nImages) ;
   imageArray87       = zeros(n_x,n_y,n_z,nImages) ;   
   imageArray85(:,:,:,1)= ps2_85_data;
   imageArray87(:,:,:,1)= ps2_87_data;
   if pars.ExpandOn == true
       imageArrayExp85    = zeros(n_x,n_r,nImages) ;
       imageArrayExp87    = zeros(n_x,n_r,nImages) ;
       imageArrayExp85(:,:,1)= ps2_85_exp;
       imageArrayExp87(:,:,1)= ps2_87_exp;
   end
   sliceCounter     = 1;
end

if saveImages2D == true && strcmp(propMode,'init') == false
    nImages            = numel(sampleTimes) + 1;
    imageArray85_xy     = zeros(n_x,n_y,nImages);
    imageArray87_xy     = zeros(n_x,n_y,nImages);
    imageArray85_zx     = zeros(n_z,n_x,nImages);   
    imageArray87_zx     = zeros(n_z,n_x,nImages);
    ps2_85_proj_xy      = squeeze(sum(ps2_85,3));
    ps2_87_proj_xy      = squeeze(sum(ps2_87,3));
    ps2_85_proj_zx      = squeeze(sum(ps2_85,1))';
    ps2_87_proj_zx      = squeeze(sum(ps2_87,1))';
    imageArray85_xy     = ps2_85_proj_xy;
    imageArray87_xy     = ps2_87_proj_xy;
    imageArray85_zx     = ps2_85_proj_zx;
    imageArray87_zx     = ps2_87_proj_zx;
    
    sliceCounter     = 1;
end   

%% CREATE FIGURES

if pars.figuresOn == true
    
% subtightplot options
stp         = [0.07,0.08];
marg_h      = 1*[0.05,0.05];
marg_w      = [0.14,0.0];

    %% Particle count
    subtightplot(4,4,8,stp);
    hold on
    plot_Nlist_85 = plot(time_vec,Nlist_85,'LineWidth',2);
    plot_Nlist_87 = plot(time_vec,Nlist_87,'LineWidth',2);
    xlim([0,Ttotal])
    grid on
    box on
    title('\textbf{Particle Count}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\textbf{Number}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    if pars.ExpandOn == true
        plot_Nlist_exp_85 = plot(time_vec,Nlist_exp_85,'LineWidth',2);
        plot_Nlist_exp_87 = plot(time_vec,Nlist_exp_87,'LineWidth',2);
        legend([plot_Nlist_85,plot_Nlist_exp_85,plot_Nlist_87,plot_Nlist_exp_87],{'Rb85','Expanded Rb85','Rb87','Expanded Rb87'})
    else
        legend([plot_Nlist_85,plot_Nlist_87],{'Rb85','Rb87'})
    end
    ylim(1.2*[0,max([N_85,N_87])])

    %% Energy
    subtightplot(4,4,4,stp);
    plot_Elist_85 = plot(time_vec,Elist_85,'LineWidth',2);
    %set(gca,'yscale','log');
    hold on
    plot_Elist_87 = plot(time_vec,Elist_87,'LineWidth',2);
    plot_Etotal   = plot(time_vec,Elist_87+Elist_85,'LineWidth',2);
    xlim([0,Ttotal])
    grid on
    title('\textbf{Mean-field Energy}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\textbf{Energy}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    legend([plot_Elist_85,plot_Elist_87,plot_Etotal],{'Rb85','Rb87','Total'})
    box on

    %% Width x
    subtightplot(4,4,13,stp)
    hold on
    plot_width_x_85 = plot(time_vec,width_x_list_85,'LineWidth',2);
    plot_width_x_87 = plot(time_vec,width_x_list_87,'LineWidth',2);
    xlim([0,Ttotal])
    grid on
    box on
    title('\textbf{Width (x)}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\textbf{Width}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    if pars.ExpandOn == true
        plot_width_r_exp_85 = plot(time_vec,width_r_exp_list_85,'LineWidth',2);
        plot_width_r_exp_87 = plot(time_vec,width_r_exp_list_87,'LineWidth',2);
        legend([plot_width_x_85,plot_width_r_exp_85,plot_width_x_87,plot_width_r_exp_87],{'Rb85','Expanded Rb85','Rb87','Expanded Rb87'})
    else
        legend([plot_width_x_85,plot_width_x_87],{'Rb85','Rb87'})
    end

    %% Width z
    subtightplot(4,4,14,stp);
    hold on
    plot_width_z_85 = plot(time_vec,width_z_list_85,'LineWidth',2);
    plot_width_z_87 = plot(time_vec,width_z_list_87,'LineWidth',2);
    xlim([0,Ttotal])
    grid on
    box on
    title('\textbf{Width (z)}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\textbf{Width}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    if pars.ExpandOn == true
        plot_width_z_exp_85 = plot(time_vec,width_z_exp_list_85,'LineWidth',2);
        plot_width_z_exp_87 = plot(time_vec,width_z_exp_list_87,'LineWidth',2);
        legend([plot_width_z_85,plot_width_z_exp_85,plot_width_z_87,plot_width_z_exp_87],{'Rb85','Expanded Rb85','Rb87','Expanded Rb87'})
    else
        legend([plot_width_z_85,plot_width_z_87],{'Rb85','Rb87'})
    end

    %% Centre-of-Mass in z
    subtightplot(4,4,15,stp);
    hold on
    plot_com_z_85      = plot(time_vec,com_z_list_85,'LineWidth',2);
    plot_com_z_87      = plot(time_vec,com_z_list_87,'LineWidth',2);
    xlim([0,Ttotal])
    grid on
    title('\textbf{Centre of Mass (z)}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
    xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    ylabel('\textbf{Position ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    if pars.ExpandOn == true
        plot_com_z_exp_85  = plot(time_vec,com_z_exp_list_85,'LineWidth',2);
        plot_com_z_exp_87  = plot(time_vec,com_z_exp_list_87,'LineWidth',2);
        legend([plot_com_z_85,plot_com_z_exp_85,plot_com_z_87,plot_com_z_exp_87],{'Rb85','Expanded Rb85','Rb87','Expanded Rb87'})
    else
        legend([plot_com_z_85,plot_com_z_87],{'Rb85','Rb87'})
    end
    box on

    %% Create unexpanded 85, 87 z-x projection
    ps2_85_proj_zx  = (abs(squeeze(sum(ps2_85,1))).^2);
    ps2_87_proj_zx  = (abs(squeeze(sum(ps2_87,1))).^2);
    
    % 85
    subtightplot(4,4,[1,2],stp);
    if useLogScale == true
        plot_density_85_zx = imagesc(x*Length*10^6,z*Length*10^6, log(1+ps2_85_proj_zx));
    else
        plot_density_85_zx = imagesc(x*Length*10^6,z*Length*10^6, ps2_85_proj_zx);
    end
    set(gca,'YDir','normal')
    ylabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    axis tight;
    colormap(jet(2048))
    plot_density_title_85_zx = title('\textbf{\boldmath$^{85}$ RbDensity, \boldmath$t = 0$}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex'); 
    
    % 87
    subtightplot(4,4,[5,6],stp);
    if useLogScale == true
        plot_density_87_zx = imagesc(x*Length*10^6,z*Length*10^6, log(1+ps2_87_proj_zx));
    else
        plot_density_87_zx = imagesc(x*Length*10^6,z*Length*10^6, ps2_87_proj_zx);
    end
    set(gca,'YDir','normal')
    ylabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    axis tight;
    colormap(jet(2048))
    plot_density_title_87_zx = title('\textbf{\boldmath$^{87}$Rb Density, \boldmath$t = 0$}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex'); 
    
    % 85, 87
    subtightplot(4,4,[9,10],stp);
    ylabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    axis tight;
    onesArray           = ones(size(ps2_85_proj_zx));
    colorArray          = ones([size(ps2_85_proj_zx),3]);
    redSlice            = log(1+ps2_85_proj_zx);
    redInd              = max(redSlice(:));
    redSlice            = redSlice/redInd;
    blueSlice           = log(1+ps2_87_proj_zx);
    blueInd             = max(blueSlice(:));
    blueSlice           = blueSlice/blueInd;
    colorArray(:,:,2)   = onesArray-redSlice/sqrt(2);
    colorArray(:,:,3)   = onesArray-redSlice/sqrt(2);
    colorArray(:,:,1)   = onesArray-blueSlice/sqrt(2);
    colorArray(:,:,2)   = colorArray(:,:,2)-blueSlice/sqrt(2);
    plot_binary_density_zx = image(x*Length*10^6,z*Length*10^6,colorArray,'CDataMapping','direct');
    set(gca,'YDir','normal')
    ylabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    plot_binary_density_title_zx = title('\textbf{Binary Density}','FontWeight','Bold','FontSize',font_LrgSize ,'interpreter','latex');
    axis tight;
    grid on
        
    %% Create unexpanded 85, 87 x,y projection
    ps2_85_proj_xy  = abs(squeeze(sum(ps2_85,3))).^2;
    ps2_87_proj_xy  = abs(squeeze(sum(ps2_87,3))).^2;
    
    % 85
    subtightplot(4,4,3,stp);
    if useLogScale == true
        plot_density_85_xy = imagesc(x*Length*10^6,y*Length*10^6, log(1+ps2_85_proj_xy));
    else
        plot_density_85_xy = imagesc(x*Length*10^6,y*Length*10^6, ps2_85_proj_xy);
    end
    set(gca,'YDir','normal')
    ylabel('\textbf{y ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    axis tight;
    colormap(jet(2048))
    plot_density_title_85_xy = title('\textbf{\boldmath$^{85}$Rb Density, \boldmath$t = 0$}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex'); 
    
    % 87
    subtightplot(4,4,7,stp);
    if useLogScale == true
        plot_density_87_xy = imagesc(x*Length*10^6,y*Length*10^6, log(1+ps2_87_proj_xy));
    else
        plot_density_87_xy = imagesc(x*Length*10^6,y*Length*10^6, ps2_87_proj_xy);
    end
    set(gca,'YDir','normal')
    ylabel('\textbf{y ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    axis tight;
    colormap(jet(2048))
    plot_density_title_87_xy = title('\textbf{\boldmath$^{87}$Rb Density, \boldmath$t = 0$}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex'); 
    
    % 85, 87
    subtightplot(4,4,11,stp);
    
    ylabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    axis tight;
    onesArray           = ones(size(ps2_85_proj_xy));
    colorArray          = ones([size(ps2_85_proj_xy),3]);
    redSlice            = log(1+ps2_85_proj_xy);
    redInd              = max(redSlice(:));
    redSlice            = redSlice/redInd;
    blueSlice           = log(1+ps2_87_proj_xy);
    blueInd             = max(blueSlice(:));
    blueSlice           = blueSlice/blueInd;
    colorArray(:,:,2)   = onesArray-redSlice/sqrt(2);
    colorArray(:,:,3)   = onesArray-redSlice/sqrt(2);
    colorArray(:,:,1)   = onesArray-blueSlice/sqrt(2);
    colorArray(:,:,2)   = colorArray(:,:,2)-blueSlice/sqrt(2);
    plot_binary_density_xy = image(x*Length*10^6,y*Length*10^6,colorArray,'CDataMapping','direct');
    set(gca,'YDir','normal')
    ylabel('\textbf{y ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    xlabel('\textbf{x ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
    plot_binary_density_title = title('\textbf{Binary Density}','FontWeight','Bold','FontSize',font_LrgSize ,'interpreter','latex');
    axis tight;
    grid on
    
    %% Create expanded z-r projection
    if pars.ExpandOn == true
        subplot(2,4,5)
        if useLogScale == true
            plot_density_exp_85 = imagesc([z,z]*Length*10^6,[-fliplr(r), r]*Length*10^6, log(1+[fliplr(abs(ps2_85_exp).^2), abs(ps2_85_exp).^2]'));
        else
            plot_density_exp_85 = imagesc([z,z]*Length*10^6,[-fliplr(r), r]*Length*10^6, [fliplr(abs(ps2_85_exp).^2), abs(ps2_85_exp).^2]');
        end
        set(gca,'YDir','normal')
        ylabel('\textbf{r ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        xlabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        axis tight;
        colormap(jet(2048))
        plot_density_exp_title_85 = title('Expanded Rb85 density @ t = 0','FontWeight','Bold','FontSize',font_LrgSize);

        subplot(2,4,6)
        if useLogScale == true
            plot_density_exp_87 = imagesc([z,z]*Length*10^6,[-fliplr(r), r]*Length*10^6, log(1+[fliplr(abs(ps2_87_exp).^2), abs(ps2_87_exp).^2]'));
        else
            plot_density_exp_87 = imagesc([z,z]*Length*10^6,[-fliplr(r), r]*Length*10^6, [fliplr(abs(ps2_87_exp).^2), abs(ps2_87_exp).^2]');
        end
        set(gca,'YDir','normal')
        ylabel('\textbf{r ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        xlabel('\textbf{z ($\mu$m)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        axis tight;
        colormap(jet(2048))
        plot_density_exp_title_87 = title('Expanded Rb87 density @ t = 0','FontWeight','Bold','FontSize',font_LrgSize);
    end



    %% Plot 2-Body interaction strength
    if RampOn == true
        subtightplot(4,4,16,stp)
        Ua_list = UaUnit85/Length^3/Energy;
        h_ua    = plot(0,Ua_list,'LineWidth',2);
        title('2-Body Interaction Strength')
        xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        ylabel('\textbf{Ua}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        xlim([0,Ttotal])
        %ylim(sort([0.9*min(rampUa_vec),1.1*max(rampUa_vec)]))
        grid on
    end
    
    %% Plot dt
    if strcmp(propMode,'prop') == true
        subtightplot(4,4,12,stp);
        hold on
        plot_dtList = plot(dtTimeList,dtList,'LineWidth',2);
        xlim([0,Ttotal])
        grid on
        title('\textbf{dt}','FontWeight','Bold','FontSize',font_LrgSize,'Interpreter','Latex');
        xlabel('\textbf{Time}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        ylabel('\textbf{Step size (s)}','FontWeight','Bold','FontSize',font_SmlSize,'Interpreter','Latex');
        set(gca,'yscale','log')
        box on 
    end

drawnow
end

%% PUSH ARRAYS TO GPU
if CUDA_flag == true
    x           = gpuArray(x);
    y           = gpuArray(y);
    z           = gpuArray(z);
    ps2_85      = gpuArray(ps2_85);
    ps2_87      = gpuArray(ps2_87);
    disp_op85   = gpuArray(disp_op85);
    disp_op87   = gpuArray(disp_op87);
    potential85 = gpuArray(potential85);
    potential87 = gpuArray(potential87);
    if bdd_on == true
        losses      = gpuArray(losses);    
    end
    if pars.ExpandOn == true
        disp_op_PostExp     = gpuArray(disp_op_PostExp);
    end
end

%% PRE-PROPAGATION
% Pre-fft and pre-Hankel state
ps2_85     = fftn(ps2_85);
ps2_87     = fftn(ps2_87);

% Begin propagation and start timer
tic

%% PROPAGATION: RAMP
fprintf('Ramping progress: ');
str1 = sprintf('00.00');
fprintf([str1,'%%']);
plotIter = 1;

% Create video
if writeVideo_flag == true && pars.figuresOn == true
    set(gcf,'color','w')
    vidObj  = VideoWriter([save_file,'_video.mj2'],'Archival');
    open(vidObj)
end
tNow    = 0;
iter = 1;
% INIT evolution
if strcmp(propMode,'init') == true
    while tNow <= Tmax + dt
        
        % Do 1 step
        propState_RK4
        
        % Update tNow
        tNow        = tNow + abs(dt);
        
        % Update plot
        if mod(plotIter,data_step) == 0
            updatePlots
            
            % Write to video
            if writeVideo_flag == true && pars.figuresOn == true
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

% PROP evolution
sampleTimeNextInd   = 1;

if strcmp(propMode,'prop') == true
    while sampleTimeNextInd <= numel(sampleTimes)
        % Get next sample time
        tNext   = sampleTimes(sampleTimeNextInd)/Time;
        sampleTimeNextInd   = sampleTimeNextInd + 1;
        getSampleFlag       = false;
        
        % Evolve system until sample time
        while getSampleFlag == false
            if RampOn == true
                Ua85    = rampUa_fun(tNow);
            end
            
            % Propagate state by 1dt
            propState_ARK45
            fprintf('Time   = %4.3e\n',tNow*Time)
            fprintf('T2Step = %4.3e\n',(tNext-tNow)*Time)
            pause(0.01)
            % Update tNow
            tNow    = tNow + abs(dt);
            
            % Update progress
%             fprintf(repmat('\b',1,length(str1)+1));
%             str1 = sprintf('%2.2f',tNow/n_t*100);
%             fprintf([str1,'%%']);
%             plotIter = plotIter+1;
        end
        
        % Get sample
        updatePlots
        
        % Write to video
        if writeVideo_flag == true && pars.figuresOn == true
            frame = getframe(main_fig);
            writeVideo(vidObj,frame);
        end
    end
end

% Close Video
if writeVideo_flag == true
    close(vidObj);
end

% UnHankel and UnFFT
ps2_85     = ifftn(ps2_85);
ps2_87     = ifftn(ps2_87);

% Renormalise
if strcmp(propMode,'init') == true
    [n_85,n_87]  = getN;
    ps2_85     = ps2_85 * sqrt(N_85/n_85);
    ps2_87     = ps2_87 * sqrt(N_87/n_87);
end

% Gather
if CUDA_flag == true
    ps2_85  = gather(ps2_85);
    ps2_87  = gather(ps2_87);
end

% it's done
time_taken = toc;
fprintf(repmat('\b',1,length(str1)+1))
fprintf('Done - ')
fprintf('Time taken: %4.3f\n',time_taken);

%% SAVE FILES
fprintf('Saving . . . ')

    %save data
    save(save_file,...
        'time_vec',...
        'width_x_list_85',...
        'width_y_list_85',...
        'width_z_list_85',...
        'com_x_list_85',...
        'com_y_list_85',...
        'com_z_list_85',...
        'width_x_list_87',...
        'width_y_list_87',...
        'width_z_list_87',...
        'com_x_list_87',...
        'com_y_list_87',...
        'com_z_list_87',...
        'Nlist_85',...
        'Nlist_87',...
        'Elist_85',...
        'Elist_87',...
        'ps2_85',...
        'ps2_87',...
        'pars',...
        'x',...
        'y',...
        'z',...
        'Length',...
        'Energy',...
        'Time');
    
    if pars.ExpandOn == true
        save(save_file,...
            'width_z_exp_list_85',...
            'width_r_exp_list_85',...
            'com_z_exp_list_85',...
            'width_z_exp_list_87',...
            'width_r_exp_list_87',...
            'com_z_exp_list_87',...
            'Nlist_exp_85',...
            '-append')
    end
    
    if saveImages3D == true && strcmp(propMode,'init') == false
        save(save_file,...
            'imageArray85',...
            'imageArray87',...
            '-append')
        if pars.ExpandOn == true
        save(save_file,...
                   'imageArrayExp85',... 
                   'imageArrayExp87',... 
                   '-append')
        end
    end
    
    if saveImages2D == true && strcmp(propMode,'init') == false
        save(save_file,...
            'imageArray85_xy',...
            'imageArray87_xy',...
            'imageArray85_zx',...
            'imageArray87_zx',...
            '-append')
    end
    
fprintf('Done\n\n')

%% AUXILIARY FUNCTIONS
    % Get potentials
    function  [p85,p87,trap85,trap87] =  getPotential(omega_x,omega_y,omega_z)
        if wgTilt==0
            zPart85     = (omega_z(1)*Time*Z).^2;
            zPart87     = (omega_z(2)*Time*Z).^2;
        else
            z_pot85     = (omega_z(1)*Time*z).^2;
            z_rot       = z*cos(wgTilt) - z_pot85*sin(wgTilt);
            z_pot_rot   = z*sin(wgTilt) + z_pot85*cos(wgTilt);
            zPart85     = interp1(z_rot,z_pot_rot,z,'spline');
            [~,~,zPart85] = meshgrid(x,y,zPart85);
            
            z_pot87     = (omega_z(2)*Time*z).^2;
            z_rot       = z*cos(wgTilt) - z_pot87*sin(wgTilt);
            z_pot_rot   = z*sin(wgTilt) + z_pot87*cos(wgTilt);
            zPart87       = interp1(z_rot,z_pot_rot,z,'spline');
            [~,~,zPart87] = meshgrid(x,y,zPart87);
            
        end
        if gravityOn == true
            yShift85        = gDQS/(omega_y(1)^2*Length);
            yShift87        = gDQS/(omega_y(2)^2*Length);
            trapCorrection  = m85*gDQS*Y*Length/Energy;
        else
            yShift85        = 0;
            yShift87        = 0;
            trapCorrection  = 0;
        end
        p85     =           0.5*((omega_x(1)*Time*X).^2 + (omega_y(1)*Time*(Y-yShift85)).^2 + zPart85)+ trapCorrection;
        p87     = m87/m85*( 0.5*((omega_x(2)*Time*X).^2 + (omega_y(2)*Time*(Y-yShift87)).^2 + zPart87)+ trapCorrection);
         
        % Use unshifted potential for ansatz
        trap85  =          0.5*((omega_x(1)*Time*X).^2 + (omega_y(1)*Time*Y).^2 + zPart85);
        trap87  = m87/m85*(0.5*((omega_x(1)*Time*X).^2 + (omega_y(1)*Time*Y).^2 + zPart87));
    end


    % get number of particles by integrating intensity
    function [n_85,n_87] = getN
        n_85 = trapz(x,trapz(y,trapz(z,abs(ps2_85).^2,3),1),2);
        n_87 = trapz(x,trapz(y,trapz(z,abs(ps2_87).^2,3),1),2);
%         if CUDA_flag == true
%            n_85     = gather(n_85);
%            n_87     = gather(n_87);
%         end
    end
    
    function [n_85,n_87] = getN_data
        n_85 = trapz(x,trapz(y,trapz(z,abs(ps2_85_data).^2,3),1),2);
        n_87 = trapz(x,trapz(y,trapz(z,abs(ps2_87_data).^2,3),1),2);
%         if CUDA_flag == true
%            n_85     = gather(n_85);
%            n_87     = gather(n_87);
%         end
    end

    % Get number of particles in expanded cloud
    function [n_85,n_87] = getExpandedN
        n_85 = trapz(z,2*pi*trapz(r,R.*abs(ps2_85_exp).^2,2),1);
        n_87 = trapz(z,2*pi*trapz(r,R.*abs(ps2_87_exp).^2,2),1);
%         if CUDA_flag == true
%            n_85     = gather(n_85);
%            n_87     = gather(n_87);
%         end
    end

    function [E_85,E_87] = getE
        % Get per-particle mean field (non-kinetic energy). Check interaction energy!!!
        warning('off','all')
        [fx,fy,fz] = gradient(ps2_85_data,x,y,z);
        KE         = 1/2*(abs(fx).^2 + abs(fy).^2 + abs(fz).^2);
        eDens_85=   real(potential85).*abs(ps2_85_data).^2 + ...
                    Ua85/2*abs(ps2_85_data).^4 + ...
                    Ua8587/2*abs(ps2_85_data).^2.*abs(ps2_87_data).^2 + ...
                    KE;
                
        %eDens_85   = KE;
                
        [fx,fy,fz] = gradient(ps2_87_data,x,y,z);
        KE         = m85/(2*m87)*(abs(fx).^2 + abs(fy).^2 + abs(fz).^2);
        eDens_87=   m87/m85*real(potential87).*abs(ps2_87_data).^2 + ...
                    Ua87/2*abs(ps2_87_data).^4 + ...
                    Ua8587/2*abs(ps2_85_data).^2.*abs(ps2_87_data).^2 + ...
                    KE;
        %eDens_87   = KE;
        %hbar^2/(2*m)*abs(gradient(ps2_data,r,z)).^2 ...   % Kinetic term. Check grad is correct
        E_85    = trapz(x,trapz(y,trapz(z,eDens_85,3),1),2);             
        E_87    = trapz(x,trapz(y,trapz(z,eDens_87,3),1),2);
    end

    function [widths,coms] = getMoments
        % Get com
        com_x_85    = trapz(x,(x).*trapz(z,trapz(y,abs(ps2_85_data).^2,1),3),2)/n85*Length*10^6;
        com_y_85    = trapz(y,(y)'.*trapz(x,trapz(z,abs(ps2_85_data).^2,3),2),1)/n85*Length*10^6;
        z_temp      = reshape(z,[1,1,numel(z)]);
        com_z_85    = trapz(z,z_temp.*trapz(y,trapz(x,abs(ps2_85_data).^2,2),1),3)/n85*Length*10^6;
        
        com_x_87    = trapz(x,(x).*trapz(z,trapz(y,abs(ps2_87_data).^2,1),3),2)/n85*Length*10^6;
        com_y_87    = trapz(y,(y)'.*trapz(x,trapz(z,abs(ps2_87_data).^2,3),2),1)/n85*Length*10^6;
        z_temp      = reshape(z,[1,1,numel(z)]);
        com_z_87    = trapz(z,z_temp.*trapz(y,trapz(x,abs(ps2_87_data).^2,2),1),3)/n85*Length*10^6;
        
        % Get widths
        width_x_85  = (1/n85)*trapz(x,(x.^2).*trapz(y,trapz(z,abs(ps2_85_data).^2,3),1),2)*Length*10^6;
        width_y_85  = (1/n85)*trapz(y,(y.^2)'.*trapz(z,trapz(x,abs(ps2_85_data).^2,2),3),1)*Length*10^6;
        z_temp      = reshape(z,[1,1,numel(z)]);
        width_z_85  = (1/n85)*trapz(z,(z_temp.^2).*trapz(x,trapz(y,abs(ps2_85_data).^2,1),2),3)*Length*10^6;
        
        width_x_87  = (1/n87)*trapz(x,(x.^2).*trapz(y,trapz(z,abs(ps2_87_data).^2,3),1),2)*Length*10^6;
        width_y_87  = (1/n87)*trapz(y,(y.^2)'.*trapz(z,trapz(x,abs(ps2_87_data).^2,2),3),1)*Length*10^6;
        z_temp      = reshape(z,[1,1,numel(z)]);
        width_z_87  = (1/n87)*trapz(z,(z_temp.^2).*trapz(x,trapz(y,abs(ps2_87_data).^2,1),2),3)*Length*10^6;
        
        % Store
        widths      = [ width_x_85,width_x_87;
                        width_y_85,width_y_87;
                        width_z_85,width_z_87];
        coms        = [ com_x_85,com_x_87;
                        com_y_85,com_y_87;
                        com_z_85,com_z_87];
                 
%         if CUDA_flag == true
%            widths   = gather(widths);
%            coms     = gather(coms);
%         end
    end
    
    % get expanded cloud width and COM 85
    function [width_z_exp, width_r_exp,com_z_exp] = getExpandedWidths85
        com_z_exp    = trapz(z,((z)'.*(2*pi*trapz(r,R.*abs(ps2_85_exp).^2,2)))/n85)*Length*10^6;
        width_z_exp  = trapz(z,((z.^2)'.*(2*pi*trapz(r,R.*abs(ps2_85_exp).^2,2)))/n85)*Length*10^6-com_z_exp^2;
        width_r_exp  = trapz([-fliplr(r) r],[-fliplr(r) r].^2.*(trapz(z,[fliplr(abs(ps2_85_exp).^2) abs(ps2_85_exp).^2]))/n85)*(Length*10^6)^2;
    end

    % get expanded cloud width and COM 87
    function [width_z_exp, width_r_exp,com_z_exp] = getExpandedWidths87
        com_z_exp    = trapz(z,((z)'.*(2*pi*trapz(r,R.*abs(ps2_87_exp).^2,2)))/n87)*Length*10^6;
        width_z_exp  = trapz(z,((z.^2)'.*(2*pi*trapz(r,R.*abs(ps2_87_exp).^2,2)))/n87)*Length*10^6-com_z_exp^2;
        width_r_exp  = trapz([-fliplr(r) r],[-fliplr(r) r].^2.*(trapz(z,[fliplr(abs(ps2_87_exp).^2) abs(ps2_87_exp).^2]))/n87)*(Length*10^6)^2;
    end

    function [g85,g87] = gpeFun(p85,p87)
        d85     = abs(p85).^2;
        d87     = abs(p87).^2;
        g85     = (1i*((-Ua85*d85-Ua8587*d87-G3*d85.^2-potential85).*p85))*dt;
        g87     = (1i*((-Ua87*d87-Ua8587*d85-G3*d87.^2-potential87).*p87))*dt;
    end

    function updateDt(newDt)
        % Update dt and dispersion operators
        dt              = newDt;
        disp_op85       = exp(-1i*0.25*dt*(Kx.^2+Ky.^2+Kz.^2));
        disp_op87       = exp(-1i*0.25*dt*(Kx.^2+Ky.^2+Kz.^2)*m85/m87);
        
        % Update dtList and dtTimeList
        dtList(end+1)       = newDt;
        dtTimeList(end+1)   = tNow*Time;
    end
    
%% PROPAGATORS
    % Propagate ps2 (combined) state by 1dt
    function propState_RK4
        % 1st half of dispersion
            % Apply operator
            ps2_85     = disp_op85.*ps2_85;
            ps2_87     = disp_op87.*ps2_87;

            % UnHankel and UnFFT
            ps2_85     = ifftn(ps2_85);
            ps2_87     = ifftn(ps2_87);

        % Interaction (RK4)
        dens85    = abs(ps2_85).^2;
        dens87    = abs(ps2_87).^2;
        A2_85     = (1i*((-Ua85*dens85-Ua8587*dens87-G3*dens85.^2-potential85).*ps2_85))*dt;
        A2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-G3*dens87.^2-potential87).*ps2_87))*dt;
        ps2_85_out= A2_85;
        ps2_87_out= A2_87;

        dens85    = abs(ps2_85+1/2*A2_85).^2;
        dens87    = abs(ps2_87+1/2*A2_87).^2;
        A2_85     = (1i*((-Ua85*dens85-Ua8587*dens87-G3*dens85.^2-potential85).*(ps2_85+1/2*A2_85)))*dt;
        A2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-G3*dens87.^2-potential87).*(ps2_87+1/2*A2_87)))*dt;
        ps2_85_out= ps2_85_out + 2*A2_85;
        ps2_87_out= ps2_87_out + 2*A2_87;

        dens85    = abs(ps2_85+1/2*A2_85).^2;
        dens87    = abs(ps2_87+1/2*A2_87).^2;
        A2_85     = (1i*((-Ua85*dens85-Ua8587*dens87-G3*dens85.^2-potential85).*(ps2_85+1/2*A2_85)))*dt;
        A2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-G3*dens87.^2-potential87).*(ps2_87+1/2*A2_87)))*dt;
        ps2_85_out= ps2_85_out + 2*A2_85;
        ps2_87_out= ps2_87_out + 2*A2_87;

        dens85    = abs(ps2_85+A2_85).^2;
        dens87    = abs(ps2_87+A2_87).^2;
        A2_85     = (1i*((-Ua85*dens85-Ua8587*dens87-G3*dens85.^2-potential85).*(ps2_85+A2_85)))*dt;
        A2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-G3*dens87.^2-potential87).*(ps2_87+A2_87)))*dt;

        ps2_85    = ps2_85+1/6*(ps2_85_out + A2_85);
        ps2_87    = ps2_87+1/6*(ps2_87_out + A2_87);
        
        % Renormalise
        if strcmp(propMode,'init') == true
            [n_85,n_87]  = getN;
            ps2_85     = ps2_85 * sqrt(N_85/n_85);
            ps2_87     = ps2_87 * sqrt(N_87/n_87);
        end
        
        % Apply lossy boundary
        if bdd_on == true
            ps2_85     = ps2_85.*losses;
            ps2_87     = ps2_87.*losses;
        end
        
        % 2nd half of dispersion
        % FFT and Hankel
        ps2_85     = fftn(ps2_85);
        ps2_87     = fftn(ps2_87);
        
        % Apply operator
        ps2_85     = disp_op85.*ps2_85;
        ps2_87     = disp_op87.*ps2_87;
    end
    
    function propState_ARK45
        ps2_85_old = ps2_85;
        ps2_87_old = ps2_87;
        
        % CHECK WHILE LOOP AND DT UPDATE WORKS
        try_again_flag = true;
        dt_changed     = false;
        while try_again_flag == true
            % 1st half of dispersion
            % Apply operator
            ps2_85     = disp_op85.*ps2_85;
            ps2_87     = disp_op87.*ps2_87;

            % UnHankel and UnFFT
            ps2_85     = ifftn(ps2_85);
            ps2_87     = ifftn(ps2_87);  

            % Interaction (ARK45)
            [A85,A87]   = gpeFun(ps2_85                                                                                                                 , ps2_87);
            [B85,B87]   = gpeFun(ps2_85 +          1/5*A85                                                                                              , ps2_87 +       1/5*A87);
            [C85,C87]   = gpeFun(ps2_85 +         3/40*A85 +       9/40*B85                                                                             , ps2_87 +       3/40*A87 +       9/40*B87);
            [D85,D87]   = gpeFun(ps2_85 +        44/45*A85 -      56/15*B85 +       32/9*C85                                                            , ps2_87 +      44/45*A87 -      56/15*B87 +       32/9*C87);
            [E85,E87]   = gpeFun(ps2_85 +   19372/6561*A85 - 25360/2187*B85 + 64448/6561*C85 - 212/729*D85                                              , ps2_87 + 19372/6561*A87 - 25360/2187*B87 + 64448/6561*C87 -   212/729*D87);
            [F85,F87]   = gpeFun(ps2_85 +    9017/3168*A85 -     355/33*B85 + 46732/5247*C85 +	49/176*D85 -    5103/18656*E85                          , ps2_87 +  9017/3168*A87 -     355/33*B87 + 46732/5247*C87 +	49/176*D87 -     5103/18656*E87);
            % 5th, 4th order estimates
            [G85,G87]   = gpeFun(ps2_85 +       35/384*A85 +          0*B85 +   500/1113*C85 + 125/192*D85 -     2187/6784*E85 +    11/84*F85           , ps2_87 +     35/384*A87 +          0*B87 +   500/1113*C87 +   125/192*D87 -     2187/6784*E87 +    11/84*F87);
            [H85,H87]   = gpeFun(ps2_85 +   5179/57600*A85 +          0*B85 + 7571/16695*C85 + 393/640*D85 -  92097/339200*E85 + 187/2100*F85 + 1/40*G85, ps2_87 + 5179/57600*A87 +          0*B87 + 7571/16695*C87 +   393/640*D87 -  92097/339200*E87 + 187/2100*F87 + 1/40*G87);

            err85 = sum(abs(G85(:)-H85(:)));%/numel(G85);
            err87 = sum(abs(G87(:)-H87(:)));%/numel(G87);
            s     = min((dt*pars.tol./(2*(Tmax)*[err85,err87])).^(1/4));
            
            % FIX PS2 ASSIGNMENTS
            if s < 1
                % Previous step was too large. Halve and try again.
                updateDt(0.5*dt)
                ps2_85  = ps2_85_old;
                ps2_87  = ps2_87_old;
                try_again_flag  = true;
                fprintf('Reducing step \n')
                
                % Reset counter
                arkStepsWithoutChange   = 0;
                
            elseif s>2 && arkStepsWithoutChange >= pars.ark_minStepsBeforeChangeDt
                % Previous step was too small. Update ps2 and double step-size for next step
                dt_changed      = true;
                if tNext-tNow < 2*dt
                    dtNew           = (tNext-tNow);
                    getSampleFlag   = true;
                    fprintf('Increasing step, but sampling \n')
                else
                    dtNew           = 2*dt;
                    fprintf('Increasing step \n')
                end
                try_again_flag  = false;
                
                
                % Reset counter
                arkStepsWithoutChange   = 0;
            else
                % Step-size was about right OR dt has been changed
                % recently
                % Update ps2 and leave dt unchanged, unless getting sample.
                if tNext-tNow < dt
                    dt_changed      = true;
                    dtNew           = (tNext-tNow);
                    getSampleFlag   = true;
                    fprintf('Keeping step, but sampling \n')
                else
                    fprintf('Keeping step \n')
                end
                try_again_flag  = false;
                
                % Increment unchanged counter
                arkStepsWithoutChange   = arkStepsWithoutChange + 1;
            end
            
            if try_again_flag == false
                % Advance
                ps2_85          = ps2_85 + G85;
                ps2_87          = ps2_87 + G87;
                
                % Renormalise
                if strcmp(propMode,'init') == true
                    [n_85,n_87]  = getN;
                    ps2_85     = ps2_85 * sqrt(N_85/n_85);
                    ps2_87     = ps2_87 * sqrt(N_87/n_87);
                end

                % Apply lossy boundary
                if bdd_on == true
                    ps2_85     = ps2_85.*losses;
                    ps2_87     = ps2_87.*losses;
                end

                % 2nd half of dispersion
                % FFT and Hankel
                ps2_85     = fftn(ps2_85);
                ps2_87     = fftn(ps2_87);

                % Apply operator
                ps2_85     = disp_op85.*ps2_85;
                ps2_87     = disp_op87.*ps2_87;
            end
            
            fprintf('dt     = %4.3e\n',dt)
        end
        % Update dt
        if dt_changed == true
            updateDt(dtNew)
        end
    end

%% FIGURE HANDLING
    % Plot data during propagations
    function updatePlots
            
            % UnHankel and UnFFT Unexpanded cloud
            ps2_85_data= ifftn(ps2_85);
            ps2_87_data= ifftn(ps2_87);
%             if CUDA_flag == true
%                 ps2_85_data= gather(ps2_85_data);
%                 ps2_87_data= gather(ps2_87_data);
%             end
            
            % Perform pre-imaging expansion
            if pars.ExpandOn == true
                % Still in Fourier/Hankel space
                ps2_85_exp     = ps2_85;
                ps2_87_exp     = ps2_85;
                
                % Perform expansion
                expandCloud                
                
                % UnHankel and UnFFT Expanded cloud
                ps2_85_exp= ifftn(ps2_85_exp);
                ps2_87_exp= ifftn(ps2_87_exp);
                if CUDA_flag == true
                    ps2_85_exp= gather(ps2_85_exp);
                    ps2_87_exp= gather(ps2_87_exp);
                end
            end
            
            % Update projections
            ps2_85_proj_zx  = gather(abs(squeeze(sum(ps2_85_data,1))).^2);
            ps2_87_proj_zx  = gather(abs(squeeze(sum(ps2_87_data,1))).^2);
            ps2_85_proj_xy  = gather(abs(squeeze(sum(ps2_85_data,3))).^2);
            ps2_87_proj_xy  = gather(abs(squeeze(sum(ps2_87_data,3))).^2);  
            
            % Store slices
            if saveImages3D == true && strcmp(propMode,'init') == false
                sliceCounter     = sliceCounter + 1;
                imageArray85(:,:,:,sliceCounter)= gather(ps2_85_data);
                imageArray87(:,:,:,sliceCounter)= gather(ps2_87_data);
                if pars.ExpandOn == true
                    imageArrayExp85(:,:,:,sliceCounter)= ps2_85_exp;
                    imageArrayExp87(:,:,:,sliceCounter)= ps2_87_exp;
                end
            end
            
            if saveImages2D == true && strcmp(propMode,'init') == false
                sliceCounter     = sliceCounter + 1;
                imageArray85_xy(:,:,sliceCounter)=ps2_85_proj_xy;
                imageArray85_zx(:,:,sliceCounter)=ps2_85_proj_zx;
                imageArray87_xy(:,:,sliceCounter)=ps2_87_proj_xy;
                imageArray87_zx(:,:,sliceCounter)=ps2_87_proj_zx;
            end
            
            % Update time
            time_vec(end+1) = tNow*Time;
            
            % Update dt
            if strcmp(propMode,'prop') == true
                set(plot_dtList,'XData',dtTimeList,'YData',dtList);
            end
            
            % Update penergy
            [E_85,E_87]        = getE;
            Elist_85(end+1)    = gather(E_85);
            Elist_87(end+1)    = gather(E_87);
%             Eplot85            = [0,diff(Elist_85)];
%             Eplot85            = sign(Eplot85).*log(abs(Eplot85));
%             Eplot87            = [0,diff(Elist_87)];
%             Eplot87            = sign(Eplot87).*log(abs(Eplot87));
%             set(plot_Elist_85,'XData',time_vec,'YData',Eplot85);
%             set(plot_Elist_87,'XData',time_vec,'YData',Eplot87);
%             set(plot_Etotal,'XData',time_vec,'YData',Eplot85+Eplot87);
            set(plot_Elist_85,'XData',time_vec,'YData',Elist_85);
            set(plot_Elist_87,'XData',time_vec,'YData',Elist_87);
            set(plot_Etotal,'XData',time_vec,'YData',Elist_87+Elist_85);
            
            % Update particle counts
            [n85,n87]          = getN_data;
            Nlist_85(end+1)    = gather(n85);
            Nlist_87(end+1)    = gather(n87);
            set(plot_Nlist_85,'XData',time_vec,'YData',Nlist_85);
            set(plot_Nlist_87,'XData',time_vec,'YData',Nlist_87);
            
            % Update Widths
                % Rb85
                [widths,coms]          = getMoments;
                width_x_list_85(end+1) = gather(widths(1,1));
                width_y_list_85(end+1) = gather(widths(2,1));
                width_z_list_85(end+1) = gather(widths(3,1));
                com_x_list_85(end+1)   = gather(coms(1,1));
                com_y_list_85(end+1)   = gather(coms(2,1));
                com_z_list_85(end+1)   = gather(coms(3,1));
                if pars.figuresOn == true
                    set(plot_width_x_85,'XData',time_vec,'YData',width_x_list_85);
                    set(plot_width_z_85,'XData',time_vec,'YData',width_z_list_85);
                    set(plot_com_z_85,'XData',time_vec,'YData',com_z_list_85);
                end

                % Rb87
                width_x_list_87(end+1) = gather(widths(1,2));
                width_y_list_87(end+1) = gather(widths(2,2));
                width_z_list_87(end+1) = gather(widths(3,2));
                com_x_list_87(end+1)   = gather(coms(1,2));
                com_y_list_87(end+1)   = gather(coms(2,2));
                com_z_list_87(end+1)   = gather(coms(3,2));
                if pars.figuresOn == true
                    set(plot_width_x_87,'XData',time_vec,'YData',width_x_list_87);
                    set(plot_width_z_87,'XData',time_vec,'YData',width_z_list_87);
                    set(plot_com_z_87,'XData',time_vec,'YData',com_z_list_87);
                end
            
            if  pars.figuresOn == true
                if useLogScale == true
                    set(plot_density_85_zx,'CData',log(1+ps2_85_proj_zx));
                    set(plot_density_87_zx,'CData',log(1+ps2_87_proj_zx));
                    set(plot_density_85_xy,'CData',log(1+ps2_85_proj_xy));
                    set(plot_density_87_xy,'CData',log(1+ps2_87_proj_xy));
                else
                    set(plot_density_85_zx,'CData',ps2_85_proj_zx);
                    set(plot_density_87_zx,'CData',ps2_87_proj_zx);
                    set(plot_density_85_xy,'CData',ps2_85_proj_xy);
                    set(plot_density_87_xy,'CData',ps2_87_proj_xy);
                end
                set(plot_density_title_85_zx,'String',sprintf('\\textbf{Rb\\boldmath$^{85}$ Density, \\boldmath$t = %4.1f$ ms}',1000*tNow*Time));
                set(plot_density_title_87_zx,'String',sprintf('\\textbf{Rb\\boldmath$^{87}$ Density, \\boldmath$t = %4.1f$ ms}',1000*tNow*Time));
                
                % Update binary mixture images
                % xz projection
                onesArray           = ones(size(ps2_85_proj_zx));
                redSlice            = log(1+ps2_85_proj_zx);
                redInd              = max(redSlice(:));
                redSlice            = redSlice/redInd;
                blueSlice           = log(1+ps2_87_proj_zx);
                blueInd             = max(blueSlice(:));
                blueSlice           = blueSlice/blueInd;
                colorArray          = zeros([size(ps2_85_proj_zx),3]);
                colorArray(:,:,2)   = onesArray-redSlice/sqrt(2);
                colorArray(:,:,3)   = onesArray-redSlice/sqrt(2);
                colorArray(:,:,1)   = onesArray-blueSlice/sqrt(2);
                colorArray(:,:,2)   = colorArray(:,:,2)-blueSlice/sqrt(2);
                set(plot_binary_density_zx,'CData',colorArray);
                
                % xy projection
                onesArray           = ones(size(ps2_85_proj_xy));
                redSlice            = log(1+ps2_85_proj_xy);
                redInd              = max(redSlice(:));
                redSlice            = redSlice/redInd;
                blueSlice           = log(1+ps2_87_proj_xy);
                blueInd             = max(blueSlice(:));
                blueSlice           = blueSlice/blueInd;
                colorArray          = zeros([size(ps2_85_proj_xy),3]);
                colorArray(:,:,2)   = onesArray-redSlice/sqrt(2);
                colorArray(:,:,3)   = onesArray-redSlice/sqrt(2);
                colorArray(:,:,1)   = onesArray-blueSlice/sqrt(2);
                colorArray(:,:,2)   = colorArray(:,:,2)-blueSlice/sqrt(2);
                set(plot_binary_density_xy,'CData',colorArray);
            end
            
            % Update Expanded Quantities
            if pars.ExpandOn == true
                % Update expanded counts
                [n1_85,n1_87] = getExpandedN;
                Nlist_exp_85(end+1)= n1_85;
                Nlist_exp_87(end+1)= n1_87;
                set(plot_Nlist_exp_85,'XData',time_vec,'YData',Nlist_exp_85);
                set(plot_Nlist_exp_87,'XData',time_vec,'YData',Nlist_exp_87);
                
                % Update Expanded Widths
                    % Rb85
                    [width_z_exp, width_r_exp,com_z_exp]  = getExpandedWidths85;
                    width_z_exp_list_85(end+1)     = width_z_exp;
                    width_r_exp_list_85(end+1)     = width_r_exp;
                    com_z_exp_list_85(end+1)       = com_z_exp;
                    set(plot_width_r_exp_85,'XData',time_vec,'YData',width_r_exp_list_85);
                    set(plot_width_z_exp_85,'XData',time_vec,'YData',width_z_exp_list_85);
                    set(plot_com_z_exp_85,'XData',time_vec,'YData',com_z_exp_list_85);
                    
                    % Rb87
                    [width_z_exp, width_r_exp,com_z_exp]  = getExpandedWidths87;
                    width_z_exp_list_87(end+1)     = width_z_exp;
                    width_r_exp_list_87(end+1)     = width_r_exp;
                    com_z_exp_list_87(end+1)       = com_z_exp;
                    set(plot_width_r_exp_87,'XData',time_vec,'YData',width_r_exp_list_87);
                    set(plot_width_z_exp_87,'XData',time_vec,'YData',width_z_exp_list_87);
                    set(plot_com_z_exp_87,'XData',time_vec,'YData',com_z_exp_list_87);
                
                % Update projections
                if useLogScale == true
                    set(plot_density_exp_85,'CData',log(1+[fliplr(abs(ps2_85_exp).^2), abs(ps2_85_exp).^2]'));
                    set(plot_density_exp_87,'CData',log(1+[fliplr(abs(ps2_87_exp).^2), abs(ps2_87_exp).^2]'));
                else
                    set(plot_density_exp_85,'CData',[fliplr(abs(ps2_85_exp).^2), abs(ps2_85_exp).^2]');
                    set(plot_density_exp_87,'CData',[fliplr(abs(ps2_87_exp).^2), abs(ps2_87_exp).^2]');
                end
                set(plot_density_exp_title_85,'String',sprintf('Expanded Rb85 Density @ t = %4.2e',tNow));
                set(plot_density_exp_title_87,'String',sprintf('Expanded Rb87 Density @ t = %4.2e',tNow));
            end
            
            % Update Ua plot
            if RampOn == true
                Ua_list(end+1)  = Ua85;
                if pars.figuresOn == true
                    set(h_ua,'XData',time_vec,'YData',Ua_list);
                end
            end
            
            if pars.figuresOn == true
                drawnow
            end
    end
    
%% CLOUD EXPANSION
    % Perform pre-imaging cloud expansion of ps2 (combined cloud)
    function expandCloud
        % Assume we're already in Fourier/Hankel space
        
        % Freespace Expansion (a_s = 0)
        strFreespace = sprintf(' - Freespace expansion progress: ');
        fprintf(strFreespace);
        strExp = sprintf('00.00');
        fprintf([strExp,'%%']);
        for expLoop = 1:n_t_exp
            
            % 1st half of dispersion
            % Apply operator
            ps2_85_exp     = disp_op.*ps2_85_exp;
            ps2_87_exp     = disp_op.*ps2_87_exp;

            % UnHankel and UnFFT
            ps2_85_exp     = iHankelFFT(ps2_85_exp);
            ps2_87_exp     = iHankelFFT(ps2_87_exp);


            % Interaction (RK4)
            dens85    = ps2_85_exp.*conj(ps2_85_exp);
            dens87    = ps2_87_exp.*conj(ps2_87_exp);
            A2_85     = (1i*((-Ua8587*dens87-G3*dens85.^2-potential).*ps2_85_exp))*dt;
            A2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-potential).*ps2_87_exp))*dt;

            dens85    = ((ps2_85_exp+1/2*A2_85).*conj(ps2_85_exp+1/2*A2_85));
            dens87    = ((ps2_87_exp+1/2*A2_87).*conj(ps2_87_exp+1/2*A2_87));
            B2_85     = (1i*((-Ua8587*dens87-G3*dens85.^2-potential).*(ps2_85_exp+1/2*A2_85)))*dt;
            B2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-potential).*(ps2_87_exp+1/2*A2_87)))*dt;

            dens85    = ((ps2_85_exp+1/2*B2_85).*conj(ps2_85_exp+1/2*B2_85));
            dens87    = ((ps2_87_exp+1/2*B2_87).*conj(ps2_87_exp+1/2*B2_87));
            C2_85     = (1i*((-Ua8587*dens87-G3*dens85.^2-potential).*(ps2_85_exp+1/2*B2_85)))*dt;
            C2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-potential).*(ps2_87_exp+1/2*B2_87)))*dt;

            dens85    = ((ps2_85_exp+C2_85).*conj(ps2_85_exp+C2_85));
            dens87    = ((ps2_87_exp+C2_87).*conj(ps2_87_exp+C2_87));
            D2_85     = (1i*((-Ua8587*dens87-G3*dens85.^2-potential).*(ps2_85_exp+C2_85)))*dt;
            D2_87     = (1i*((-Ua87*dens87-Ua8587*dens85-potential).*(ps2_87_exp+C2_87)))*dt;

            ps2_85_exp     = ps2_85_exp+1/6*(A2_85+2*B2_85+2*C2_85+D2_85);
            ps2_87_exp     = ps2_87_exp+1/6*(A2_87+2*B2_87+2*C2_87+D2_87);

            % Apply lossy boundary
            if bdd_on == true
                ps2_85_exp     = ps2_85_exp.*losses;
                ps2_87_exp     = ps2_87_exp.*losses;
            end
        
            % 2nd half of dispersion
            % FFT and Hankel
            ps2_85_exp     = HankelFFT(ps2_85_exp);
            ps2_87_exp     = HankelFFT(ps2_87_exp);

            % Apply operator
            ps2_85_exp     = disp_op.*ps2_85_exp;
            ps2_87_exp     = disp_op.*ps2_87_exp;
            
            % Update progress
            fprintf(repmat('\b',1,length(strExp)+1));
            strExp = sprintf('%2.2f',expLoop/n_t_exp*100);
            fprintf([strExp,'%%']);
        end
        fprintf(repmat('\b',1,length(strExp)+length(strFreespace)+1))
        
        % Post-Freespace Expansion (a_s != 0)
        strPostFreespace = sprintf(' - Post-freespace expansion progress: ');
        fprintf(strPostFreespace);
        strPExp = sprintf('00.00');
        fprintf([strPExp,'%%']);
        for expLoop = 1:n_t_postexp
            
            % 1st half of dispersion
            % Apply operator
            ps2_85_exp     = disp_op.*ps2_85_exp;
            ps2_87_exp     = disp_op.*ps2_87_exp;
            
            % UnHankel and UnFFT
            ps2_85_exp     = iHankelFFT(ps2_85_exp);
            ps2_87_exp     = iHankelFFT(ps2_87_exp);
            
            
            % Interaction (RK4)
            dens85    = ps2_85_exp.*conj(ps2_85_exp);
            dens87    = ps2_87_exp.*conj(ps2_87_exp);
            A2_85     = (1i*((-Ua_exp_85*dens85-Ua_exp_8587*dens87-G3*dens85.^2-potential).*ps2_85_exp))*dt;
            A2_87     = (1i*((-Ua87*dens87-Ua_exp_8587*dens85-potential).*ps2_87_exp))*dt;
            
            dens85    = ((ps2_85_exp+1/2*A2_85).*conj(ps2_85_exp+1/2*A2_85));
            dens87    = ((ps2_87_exp+1/2*A2_87).*conj(ps2_87_exp+1/2*A2_87));
            B2_85     = (1i*((-Ua_exp_85*dens85-Ua_exp_8587*dens87-G3*dens85.^2-potential).*(ps2_85_exp+1/2*A2_85)))*dt;
            B2_87     = (1i*((-Ua87*dens87-Ua_exp_8587*dens85-potential).*(ps2_87_exp+1/2*A2_87)))*dt;
            
            dens85    = ((ps2_85_exp+1/2*B2_85).*conj(ps2_85_exp+1/2*B2_85));
            dens87    = ((ps2_87_exp+1/2*B2_87).*conj(ps2_87_exp+1/2*B2_87));
            C2_85     = (1i*((-Ua_exp_85*dens85-Ua_exp_8587*dens87-G3*dens85.^2-potential).*(ps2_85_exp+1/2*B2_85)))*dt;
            C2_87     = (1i*((-Ua87*dens87-Ua_exp_8587*dens85-potential).*(ps2_87_exp+1/2*B2_87)))*dt;
            
            dens85    = ((ps2_85_exp+C2_85).*conj(ps2_85_exp+C2_85));
            dens87    = ((ps2_87_exp+C2_87).*conj(ps2_87_exp+C2_87));
            D2_85     = (1i*((-Ua_exp_85*dens85-Ua_exp_8587*dens87-G3*dens85.^2-potential).*(ps2_85_exp+C2_85)))*dt;
            D2_87     = (1i*((-Ua87*dens87-Ua_exp_8587*dens85-potential).*(ps2_87_exp+C2_87)))*dt;
            
            ps2_85_exp     = ps2_85_exp+1/6*(A2_85+2*B2_85+2*C2_85+D2_85);
            ps2_87_exp     = ps2_87_exp+1/6*(A2_87+2*B2_87+2*C2_87+D2_87);
            
            % Apply lossy boundary
            if bdd_on == true
                ps2_85_exp     = ps2_85_exp.*losses;
                ps2_87_exp     = ps2_87_exp.*losses;
            end
            
            % 2nd half of dispersion
            % FFT and Hankel
            ps2_85_exp     = HankelFFT(ps2_85_exp);
            ps2_87_exp     = HankelFFT(ps2_87_exp);
            
            % Apply operator
            ps2_85_exp     = disp_op.*ps2_85_exp;
            ps2_87_exp     = disp_op.*ps2_87_exp;
            
            
            % Update progress
            fprintf(repmat('\b',1,length(strPExp)+1));
            strPExp = sprintf('%2.2f',expLoop/n_t_postexp*100);
            fprintf([strPExp,'%%']);
        end
        %fprintf(repmat('\b',1,length(strPExp)+1))
        fprintf(repmat('\b',1,length(strPExp)+ +length(strPostFreespace)+1))

    end

end

%% EXTERNAL FUNCTIONS
% Send variables back to main workspace
function evalstrings = assignin_base(vars)
evalstrings = cell(numel(vars),1);

for i = 1:numel(vars)
    str             = ['assignin(''base'',''',vars{i},''',',vars{i},')'];
    evalstrings{i}  = str;
end
end

% Create absorbing 'boundaries' using 40th order super-gaussian.
% Attenuates everywhere, however effect is small in ROI. Small attenuation
% suppresses resonances of high frequencies, without affecting atom number
% significantly
function mask = boundaries(n_i,n_j,n_k)
    x = linspace(-1.1,1.1,n_i);
    y = linspace(-1.1,1.1,n_j);
    z = linspace(-1.1,1.1,n_k);
    [X,Y,Z] = meshgrid(x,y,z);
    mask = exp(-(X.^40+Y.^40+Z.^40));
end
