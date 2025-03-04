%   This code implements the state-space multitaper spectrogram  
%   described in Kim et al., 2018 PNAS. 
%
%   Usage:
%   main.m: Main code    
%   EM_parameters.m: Compute noise & state variance using EM algorithm 
%   periodogram.m: Compute periodogram
%   multitaper.m: Compute multitaper spectrogram
%   SS_ST.m: Compute SS periodogram
%   SS_MT.m: Compute SS mutitaper spectrogram
%
%   From the paper:
%  "State-space multitpaer time-freqeuncy analysis"
%   Kim, S-E, Behr, MK, Ba, D & Brown, EN
%   PNAS, 2018
%
%   Copyright 2018 The General Hospital Coporation, authored by Seong-Eun Kim, Ph.D.
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Last modified 1/11/2018
%
%************************************************************************** 

function [spect, spect_taper, fourier_coefficients] = multitaper_fc(yy, fs, TW, K)

[nw, N] = size(yy); %nw=frequency bins, N=time windows, K=tapers

[tapers,concentrations]=dpss(nw,TW,K);

spect = zeros(nw,N);
spect_taper = zeros(nw,K,N);
fourier_coefficients = zeros(nw,K,N); % Preallocate memory
y_ex = repmat(yy,1,1,K);
y_ex = permute(y_ex,[1 3 2]);
for i = 1 : N
%     J = fft(tapers.*y_ex(:,:,i),nw)/sqrt(nw);
    J = fft(tapers.*y_ex(:,:,i),nw);
    fourier_coefficients(:,:,i) = J/fs; % Store complex Fourier coefficients
    spect_taper(:,:,i) = (conj(J).*J)/fs;
    spect(:,i)=squeeze(mean(conj(J).*J,2))/fs;
end