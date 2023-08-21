% To run the script download example files from 
% https://lu.app.box.com/folder/222563865593?s=9zd90x2rd2krexwayhhe7zc0c9mnxjau
% there you can find ten simulation results from a Vienna scenario for V2V communication

clear
close all
%% Article: Laura Bernado ́ "Delay and Doppler Spreads of Non-Stationary Vehicular Channels for Safety Relevant Scenarios"
% Laura Bernado's channel sounding parameters:
LB_t_sound = 307.2e-6; % snapshots periodicity in sec => maximum Doppler shift ~ 1600 Hz
LB_B_sound = 240e6;    % 240 MHz
LB_F_sound = 5.6e9;    % 5.6 GHz
LB_Q_sound = 769;      % subcarriers H(t,f)
LB_S_sound = 32500;    % 32500 snapshots
LB_f_sound = LB_B_sound/LB_Q_sound;% subcarrier spacing
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% H(m,q): time index m goes [0,...,S-1]; frequency index q goes [0,...,Q-1]
%
%//////////////////////////////////////////////////////////////////////////
% Magic numbers I and J:
% I orthogonal time-domain tapers and J orthogonal frequency-domain tapers
% resulting in a total of IJ orthogonal two-dimensional tapering functions
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% LSF parameters:
% fading process is locally stationary within MxN samples in time and freq
%
%//////////////////////////////////////////////////////////////////////////
LB_M_sound = 128;
LB_N_sound = 128;

LB_tau_sound = LB_Q_sound/(LB_B_sound*LB_N_sound);  % ~ 25 nsec
LB_Dop_sound = 1/(LB_t_sound*LB_M_sound);        % ~ 25.43 Hz
LB_avg_sound = LB_t_sound*LB_M_sound;            % reverse if Dop_sound
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Our aim in the simulation is to find parameters M and N to match with the
% sounding tau_sound and Dop_sound parameters since the speeds are the same
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Read measurements data: example data can be found in https://lu.app.box.com/folder/222563865593?s=9zd90x2rd2krexwayhhe7zc0c9mnxjau
% It is better to save different formats of files in the .mat format
filesFolder = 'ViennaScenarioUnitySimulationResults/'; % choose your own folder

% In the current version, the program process SIVERT1 simulation results from
% a Vienna Scenario described in Markus Hoffer paper https://ieeexplore.ieee.org/document/10012965


fileID = 1;


filePath = [filesFolder, 'matData', num2str(fileID), '.mat'];
load(filePath); % Output matrix HftData
% Number of snapshots
nt = size(HftData,2);
%% Simulation
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% Simulation parameters:
%
%//////////////////////////////////////////////////////////////////////////
ChaSi_F_sim = 3.2e9;      % 3.2 GHz
ChaSi_t_sim = 1600e-6;    % 1600 usec
ChaSi_Q_sim = 311;        % subcarriers H(t,f)
ChaSi_f_sim = 500e3;      % 500 kHz
ChaSi_B_sim = ChaSi_Q_sim*ChaSi_f_sim;

ChaSi_M_sim = 90;%round(LB_avg_sound/(ChaSi_t_sim)); % basically, M_sim does not mean a lot in
% terms of Nyquist theorem and we can relax the requirements for M
ChaSi_N_sim = ChaSi_Q_sim; % we have no assumption regarding stationarity properties in the
% frequency domain

% duration of sounding after averaging
ChaSi_avg_sound = ChaSi_M_sim*ChaSi_t_sim;

% maximum captured speed
c0 = physconst('LightSpeed');
v_max = c0/ChaSi_F_sim*(1/(2*ChaSi_t_sim));
disp(['Maximum speed is ', num2str(v_max),' m/s'])

% calculate how many time instances kt we have
ChaSi_N_kt = floor(nt/ChaSi_M_sim);
% calculate time steps after widowing
ChaSi_kt_time = (1:ChaSi_N_kt)*(ChaSi_M_sim*ChaSi_t_sim);
% calculate delay steps
dt = 1/(ChaSi_Q_sim*ChaSi_f_sim);
ChaSi_tau_steps = (dt:dt:ChaSi_Q_sim*dt)*1e6;
% calculate Doppler steps
ChaSi_dop_steps = linspace(-v_max,v_max,ChaSi_M_sim);


% H matrix splitting into slices with the hight ChaSi_N_sim and width ChaSi_M_sim
H_partitioned = reshape(HftData(:, 1:ChaSi_N_kt*ChaSi_M_sim), ChaSi_N_sim, ChaSi_M_sim, ChaSi_N_kt);

%% Generate DPSS matrices
I = 2;
J = 1; % from Markus Hoffer paper

% the results of the builting function dpss give results according to
% the Slepian theory considering the bandlimit (p. 33 Laura Bernado thesis)
dpssSeqI = dpss(ChaSi_M_sim, I); dpssSeqI = dpssSeqI(:,1:I);
dpssSeqJ = dpss(ChaSi_N_sim, J); dpssSeqJ = dpssSeqJ(:,1:J);
G_mtrx = zeros(ChaSi_N_sim, ChaSi_M_sim,I*J);
for i = 1:I
    for j = 1:J
        w = (i-1)*J + j;
        G_mtrx(:,:, w) = dpssSeqJ(:,j)*dpssSeqI(:,i)';
    end
end

%% LSF Estimation

% preparing DFT matrices
if mod(ChaSi_M_sim, 2) == 0
    p = -ChaSi_M_sim/2:ChaSi_M_sim/2-1;
else
    p = -(ChaSi_M_sim-1)/2:(ChaSi_M_sim-1)/2;
end
m = p';
FT_Mtrx = exp(-1j*2*pi*m*p/ChaSi_M_sim);
% in all papers, starting from 2008, they provide wrong formulas
% thus we implement a correct FT
if mod(ChaSi_N_sim, 2) == 0
    q = -ChaSi_N_sim/2:ChaSi_N_sim/2-1;
else
    q = -(ChaSi_N_sim-1)/2:(ChaSi_N_sim-1)/2;
end
n = (0:ChaSi_N_sim-1)';
IFT_Mtrx = exp(1j*2*pi*n*q/ChaSi_N_sim);
% in all papers, starting from 2008, they provide wrong formulas
% thus we implement a correct FT by deviding 2pi by the size of FT

C_LSF = zeros(ChaSi_N_sim, ChaSi_M_sim, ChaSi_N_kt);
P_tau = zeros(ChaSi_N_sim, ChaSi_N_kt);
P_dop = zeros(ChaSi_N_kt, ChaSi_M_sim);

for kt = 1:ChaSi_N_kt
    for w = 1:J*I
        tempM = H_partitioned(:,:,kt).*G_mtrx(:,:,w);
        Hktnp = IFT_Mtrx*tempM*FT_Mtrx/ChaSi_M_sim;
        C_LSF(:,:,kt) = C_LSF(:,:,kt) + (Hktnp.*conj(Hktnp))/(J*I);
    end

    P_tau(:,kt) = sum(C_LSF(:,:,kt),2)/ChaSi_M_sim;
    P_dop(kt,:) = sum(C_LSF(:,:,kt),1)/ChaSi_N_sim;
end

disp('LSF has been calculated')

%% plotting figures

% PDP
hfig1 = figure(fileID);
surf(ChaSi_kt_time, ChaSi_tau_steps, 10*log10(P_tau/max(max(P_tau))),'linestyle','none')
shading flat;
view([0 90]);
xlim([0 max(ChaSi_kt_time)]);
ylim([0 2]);
yticks([0, 0.5, 1, 1.5, 2])
caxis([-60 0]);
xsizepixels = 750; % 560
ysizepixels = 450; % 420

ylabel('delay [us]');
xlabel('time [s]');
colorbar
title({['PDP ' num2str(3.2) ' GHz']});
colormap hot
set(gca,'FontSize',16)
fig=gcf;
fig.Position(3:4)=[xsizepixels,ysizepixels];


% DSD
hfig2 = figure(2);
surf(ChaSi_kt_time, ChaSi_dop_steps,10*log10(P_dop'/max(max(P_dop))),'linestyle','none')
shading flat;
view([0 90]);
xlim([0 max(ChaSi_kt_time)]);
ylim([min(ChaSi_dop_steps), max(ChaSi_dop_steps)]);
yticks([-20 -10 0 10 20])
caxis([-60 0]);
xsizepixels = 750; % 560
ysizepixels = 450; % 420

ylabel('rel. velocity [m/s]');
xlabel('time [s]');
colorbar
title({['DSD ' num2str(3.2) ' GHz']});
colormap hot
set(gca,'FontSize',16)
fig=gcf;
fig.Position(3:4)=[xsizepixels,ysizepixels];

