function AudioPower_SeriesOfScalar = AP(auData,totalSampleNum,samplingRate,scalingRatio,elementNum,weight) 

%% File: AudioPowerType.m
%%
%% ------------- AudioPowerType--------------------
%%
%% The function of this subroutine is to describe the temporally-smoothed 
%% instantaneous power (square of waveform values).
%%
%% 
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

%% Example: AudioPower_SeriesOfScalar = AP(A,size(A),B,4,64,[]);
%%          plot(AudioPower_SeriesOfScalar);
                
if length(weight)==0
   weight_flag = 0;
else
   weight_flag = 1;
end
hopSize = samplingRate;
sampleNum = 1024;
frameNum = floor(totalSampleNum/sampleNum);

sumElement = sum(elementNum);
AudioPower_SeriesOfScalar = zeros(frameNum, sumElement);
for i = 1:frameNum
   signal = auData(1+(i-1)*sampleNum:i*sampleNum);
   audioPowerData = signal.^2;
   meanValues(i,:) = h_Mean_SeriesOfScalar(audioPowerData, scalingRatio, elementNum, weight_flag,weight, 0);
end

AudioPower_SeriesOfScalar = meanValues;