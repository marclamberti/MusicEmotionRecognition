%% GIA NA TREKSEI: 
%% [A,B,C] = wavread('chimes.wav');
%% [H,U] = AH(A,length(A)/C,B); 
%% plot(H); stem(U);

%% TI KANEI: Ypologizei ton vathmo armonikotitas enos simatos. San eksodo
%% exei to harmonicRatio, to opoio einai i perigrafi tou logou tis
%% armonikis isxyos pros tin sinoliki isxy, kai to upperLimitOfHarmonicity,
%% to opoio einai oi syxnotites pera apo autes to fasma den mporei na
%% thewrithei armoniko (1 sixnotita gia kathe frame).

function [harmonicRatio,upperLimitOfHarmonicity] = AudioHarmonicity(auData,totalSampleNum,samplingRate) 

%% File: AudioSpectrumSpread.m
%% 
%% ----------------AudioHarmonicityType-----------------------
%% The function of this subroutine is to describe 
%% the degree of harmonicity of an audio signal 
%% AudioHarmonicityType is a description of the 
%% spread of the log-frequency power spectrum.
%%
%% input:
%% -auData:                          incoming signal
%% -totalSampleNum:                  total Sample Number of signal
%% -samplingRate:                    sampling rate of the signal
%% 
%% output:
%% -harmonicRatio:                   Series of values of the harmonic ratio
%% -upperLimitOfHarmonicity:         Series of values of the UpperLimitOfHarmonicity
%% 

%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)
% === Initialization  ===

hopSize = 1;
sampleNum = fix(samplingRate*0.032);
sampleNum1 = fix(samplingRate*0.04);
frameNum = floor(totalSampleNum/sampleNum1);

% === analysis parameters  ===

windowSize = sampleNum;
window = hanning(windowSize);
fftSize = 2^nextpow2(windowSize);
freqResolution = samplingRate/fftSize;
lowFreqNum = fix(62.5/freqResolution);

upperLimitOfHarmonicity =[];    
   harmonicRatio =[];
%harmonicRatio = zeros(frameNum,fftSize/2-1);
%upperLimitOfHarmonicity = zeros(1, frameNum);
%x =zeros(1,fftSize/2-1);

for i = 1:frameNum
   signal = auData(1+(i-1)*sampleNum1:i*sampleNum1);
   %% === step a: calculate r(k), minj_pos(r(k)), h_k ===
   [minj_pos, H_k] = h_PeriodicSignalDetection(signal,samplingRate);
   if H_k <= 1.e-4
   	harmonicRatio(i,:) = zeros(1, sampleNum);
   	upperLimitOfHarmoncity(i) = 0;
   	return
   else
      signal = auData(1+(i-1)*sampleNum:i*sampleNum);
      [minj_pos, H_k] = h_PeriodicSignalDetection(signal,samplingRate);
      if i*(sampleNum+ ceil(minj_pos)) <= totalSampleNum
         signal1 = auData(1+(i-1)*(sampleNum+ceil(minj_pos)):i*(sampleNum+ceil(minj_pos)));
      signal2 = signal1';
      else 
         signal1 = [signal' zeros(1,ceil(minj_pos))];
      signal2 = signal1;
      end
      %signal2 = signal1';
      %[harmonicRatio(i,:),upperLimitOfHarmoncity(i)] = harmonicity(signal2,window,samplingRate,minj_pos);
      %% === Memory Allocation ==
		snCombfiltered = zeros(sampleNum,1);
		
		%% === Initialization ===
		harmonicityRatio = [];
		lowHarmonicity = [];
		aflim = [];
	   lowPos = [];  
      %% === Calculate HarmonicRatio ===
      %% === step b) calculate DFT of s(i) and s(i)-s(i+j) === 
		sn = signal2(1:sampleNum);
        snFFT = h_amplFFT(sn,window');
		snPower = snFFT.^2;
		for i=1:sampleNum
   		 pos1 = i+minj_pos;
   		 sniplusj = interp1(signal2, pos1); 
         if length(sniplusj)~=0 
   		    snCombfiltered(i) = signal2(i)-sniplusj;
         else 
            snCombfiltered(i) = signal2(i);
         end
        end
		snCombFiltFFt = h_amplFFT(snCombfiltered,window);
		snCombPower = snCombFiltFFt.^2;

		%% === Calculate power spectra and group the components below 62.5 Hz ===
		%%fftSize = length(snCombFiltFFt);
		%%freqResolution = samplingRate/fftSize;
		%%loFreqNum = fix(62.5/freqResolution);
		harmonicityRatio =[harmonicityRatio sum(snCombPower(1:lowFreqNum))/sum(snPower(1:lowFreqNum))];

		%% === For each frequency f calculate the sum of power beyond that frequency,===
		%% === for both the original and comb-filtered signal, and take their ratio. ===

		for k=lowFreqNum+1:fftSize/2
   		af = sum(abs(snCombPower(k:fftSize/2)).^2)/sum(abs(snPower(k:fftSize/2)).^2);
   		aflim = [aflim af];
		end
		harmonicityRatio = [harmonicityRatio aflim];

		%% === Starting from fmax(sampleNum)  and moving down in frequency, find the ===
		%% === lowest frequency for which this ratio is smaller than a threshold (0.5). ===
		%% === If that frequency is 0, replace it by 31.25 Hz   ===

		for k=length(aflim):-1:1
 			if (aflim(k) < 0.5)
   			lowHarmonicity = [lowHarmonicity aflim(k)];
   			lowPos = [lowPos (k-1)];
 			end
        end
		
        if length(lowPos)~=0 
        j = lowPos(length(lowPos));
        else j=0;
        end
		if j == 0 & length(j) == 0
    		fmin = 31.25;
		else
    fmin = j*freqResolution;
	end

	%% === Convert this value to an octave scale based on 1 kHz. ? ===
	freq1KIndex = ceil(1000/freqResolution);
	freqUpper = ceil((samplingRate/2)/freqResolution);
   upperLimitOfHarmonicity =[upperLimitOfHarmonicity fmin];    
   harmonicRatio =[harmonicRatio harmonicityRatio];

   %upperLimitOfHarmonicity(i)= x
   %harmonicRatio(i,:) =y 
%upperLimitOfHarmonicity  

   %	harmonicRatio(i,:) = x;
    %  UpperLimitOfHarmoncity(i) = y;
   end
end