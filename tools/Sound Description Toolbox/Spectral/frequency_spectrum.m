% Gia na treksei: freqspectrum = frequency_spectrum('FAC_0A.wav',0.16);

% Ti kanei: Einai mia boithitiki sinartisi.
% Kanoume framing sto shma me parauthiro hamming.
% thetoyme to win overlap sth timh winLen/10 gia logous katanalosis mnimis 
% (den eparkouse akoma kai gia sxetika mikra simata) alla kai gia logoys
% taxititas. Sti sinexeia gai kathe frame tou simatos kanoyme DFT afou
% prota frontisoume o arimos ton stoixeion kathe frame na einai dynami toy
% 2 (an den einai gemizoume me midenika).

function  freqspectrum = frequency_specrtum(filename,frameperiod) 

% Default mpeg-7 value for frame period is 0.016

[speechSignal, Fs, Bits] = wavread(filename);

% A hamming window is chosen
winLen = floor(length(speechSignal)*frameperiod);
winOverlap = floor(winLen/300);
wHamm = hamming(winLen);

% Framing and windowing the signal without for loops.
sigFramed = buffer(speechSignal, winLen, winOverlap, 'nodelay');
sigWindowed = diag(sparse(wHamm)) * sigFramed;

% frequency_specrtum calculation
[SampleNumPerFrame, TotalFrameNum]=size(sigWindowed);
freqspectrum=[];
n=2^nextpow2(SampleNumPerFrame);
for i=1:TotalFrameNum
freqspectrum=[freqspectrum  abs(fft(sigWindowed(:,i),n))];
end
