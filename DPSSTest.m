% This scripts is intended for the dpss theory implementation sanity check
% Discrete Prolate Spheroidal (Slepian) Sequences - dpss

% What it does:
% This script explains the relationship between the output of the Matlab 
% dpss function and Slepian sequences. More details about Slepian sequences
% you can find in the Slepian's original paper "Prolate spheroidal wave 
% functions, fourier analysis, and uncertainty — V: the discrete case"

% Matlab has a built-in function for calculating dpss sequences. But often,
% Matlab users cannot understand how Matlab generates these sequences and
% how these sequences are related to the Slepian sequences.


clear
close all

%% Defining DPSS sequences parameters
% Parameter M is often defined by the Fourier transform size/window.
M = 311;

% define grid
l = 0:M-1;
m = 0:M-1;

% An unclear parameter I most probably this parameter is related to the 
% bandwidth of the signal under investigation.
I = 2;

%% searching for eigen values from the Toeplitz matrix build from specific
% SINC functions defined by parameters M and I. More details in the paper
% Prolate Spheroidal Wave Functions. Fourier Analysis, and Uncertainty—V:
% The Discrete Case || by D. SLEPIAN
S1 = zeros(M);
for j = 1:M
    for i = 1:M
        if l(i) ~= m(j)
            S1(j,i) = sin( 2*pi*(I/M)*( l(i) - m(j) ) )/(pi*( l(i) - m(j) ));
        else
            S1(i,j) = 2*(I/M);
        end
    end
end
[V,D] = eig(S1);
eigD = diag(D)';
% drawing the Toeplitz matrix
surf(S1,'LineStyle','none')

%% Obtaining sequences from the Matlab built-in dpss function
% running Matlab function to compare the results with the theoretical results
[I_dpss_tapers,I_lambda] = dpss(M, I);

%% Visualization results
% drawing prolate spheroidal sequences obtained from eigen vectors of
% Toeplitz matrix
figure
plot(flipud(V(:,end-(I-1):end)),'g','LineWidth',2)
hold on
% and from the Matlab built-in dpss function
plot(I_dpss_tapers(:,1:I),'k')


