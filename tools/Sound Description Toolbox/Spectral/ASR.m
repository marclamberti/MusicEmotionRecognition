% Audio Spectrum Rolloff Frequency Estimator
% Example: AudioSpectrumRolloff = ASR('BaCl.mf.C4B4_3.wav',0.016);

% Description: Einai mia metriki toy se poso upsiles sixnotites sto fasma
% iparxei ena sigkekrimeno tmima tiw energeias.

function AudioSpectrumRolloff = ASR(filename,frameperiod) 

%change of number-format
format long;
 
%default value of the threshold TH (range 0-1)
TH=0.92;

%calling the help-function frequency_spectrum
freqspectrum = frequency_spectrum(filename,frameperiod);
[SampleNumPerFrame, TotalFrameNum]=size(freqspectrum);

% Right part initialization
tempSamplesSum = 0;
tempFramesSum= [];

% Right Part of the equation calculation
for n=1:TotalFrameNum
    for k=1:SampleNumPerFrame
        tempSamplesSum =  freqspectrum(k,n) + tempSamplesSum;
    end  
    tempFramesSum(n) = TH * tempSamplesSum;
    tempSamplesSum = 0;
end

% Left part initialization
tempSamplesSum2 = 0;
tempFramesSum2= [];

for n=1:TotalFrameNum
    AudioSpectrumRolloff(n) = SampleNumPerFrame;
end    

% Left Part of the equation calculation and checking the condition
for n=1:TotalFrameNum
    for k=1:SampleNumPerFrame
        tempSamplesSum2 =  freqspectrum(k,n) + tempSamplesSum2;
        if tempSamplesSum2 < tempFramesSum(n) 
           AudioSpectrumRolloff(n)=k; 
        end;                         
    end  
    tempSamplesSum2 = 0;
end
