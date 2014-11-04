function amplFft_SeriesOfScalar =h_amplFFT(signal,window) 

%% file: amplFFT.m 
%% 
%% The function of this file is to compute the magnitude
%% of a signal in frequency domain
%%
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

%% === analysis parameters
windowSize = length(window);
fftSize = 2^nextpow2(windowSize);

% === amplitude STFT ===
data = signal.*window;
amplFft_SeriesOfScalar = abs(fft(data,fftSize));
