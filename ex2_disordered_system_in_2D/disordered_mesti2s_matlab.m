%% Open Channel Through Disorder
% In this example, we show how to use mesti2s() to compute the scattering
% matrix of a strongly scattering disordered medium, analyze the scattering
% matrix for the TM polarization.

clear all

%% System parameters

% Dimensions of the system, in units of the wavelength lambda_0
dx      = 1/10;  % discretization grid size
W       = 30;    % width of the scattering region
L       = 12;    % thickness of the scattering region
L_tot   = 40;    % full length of the system for plotting
r_min   = 0.2;   % minimal radius of the cylindrical scatterers
r_max   = 0.4;   % maximal radius of the cylindrical scatterers
min_sep = 0.05;  % minimal separation between cylinders
nummber_density = 1.3; % number density, in units of 1/lambda_0^2
rng_seed = 0;    % random number generator seed

% Relative permittivity, unitless
epsilon_scat = 2.0^2; % cylindrical scatterers
epsilon_bg   = 1.0^2; % background in the scattering region
epsilon_L    = 1.0^2; % frees space on the left
epsilon_R    = 1.0^2; % frees space on the right

yBC = 'periodic'; % boundary condition in y

% Generate a random collection of non-overlapping cylinders
% Note subpixel smoothing is not applied for simplicity
build_TM = true;
build_TE = true;
[epsilon, inv_epsilon, x0_list, y0_list, r0_list] = ...
    build_epsilon_disorder(W, L, r_min, r_max, min_sep, nummber_density, ...
    rng_seed, dx, epsilon_scat, epsilon_bg, build_TM, build_TE, yBC);


%% Compute the transmission matrix for TM polarization
syst.polarization = 'TM';
syst.epsilon = epsilon;

syst.epsilon_L = epsilon_L;
syst.epsilon_R = epsilon_R;
syst.length_unit  = 'lambda_0';
syst.wavelength = 1;
syst.dx = dx;
syst.yBC = yBC;

PML.npixels_PML = 20;
syst.xBC={PML,PML};

[S_matlab,channels] = mesti2s(syst, {'left','right'}, {'left','right'});

epsilon_zz = 1./inv_epsilon{1};
epsilon_yy = epsilon;
epsilon_xx = 1./inv_epsilon{2};

% Save the mat files for later comparsion
save('epsilon_xx.mat','epsilon_xx')
save('epsilon_yy.mat','epsilon_yy')
save('epsilon_zz.mat','epsilon_zz')
save('S_matlab.mat','S_matlab')
