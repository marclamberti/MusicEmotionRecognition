% function [ZeroCrossingRate avgZeroCrossingRate] = ZCR(speechSignal,frameperiod) 

% Example: 
% [A,B,C] = wavread('BaCl.mf.C4B4_3.wav');
% [ZCR avZCR] = ZCR(A,0.016);

% DESCRIPTION: Computes the Zero-Crossing Rate and the average ZCR of a signal.
% Praktika to Zero Crossing Rate einai o ariumos twn midenismwn tou simatos
% sti monada tou xronou. Kanoume framing sto shma me parathiro hamming.
% thetoyme to win overlap sth timh winLen/10 gia logous katanalwsis mnimis 
% (den eparkouse akoma kai gia sxetika mikra simata) alla kai gia logoys
% taxititas. 

function [ZeroCrossingRate avgZeroCrossingRate] = ZCR(speechSignal,frameperiod) 
% Default mpeg-7 value for frame period is 0.016

% A hamming window is chosen
winLen = floor(length(speechSignal)*frameperiod);
winOverlap = floor(winLen/10);
wHamm = hamming(winLen);

% Framing and windowing the signal without for loops.
sigFramed = buffer(speechSignal, winLen, winOverlap, 'nodelay');
sigWindowed = diag(sparse(wHamm)) * sigFramed;

% Zero-Crossing Rate calculation
[SampleNumPerFrame, TotalFrameNum]=size(sigWindowed);
M=SampleNumPerFrame;
for i=1:TotalFrameNum
    for j=2:SampleNumPerFrame-1
        TempSignDifference(i,j)=(1/(SampleNumPerFrame-1)) * sum (abs(sign(sigWindowed(j,i)) - sign(sigWindowed(j-1,i))));
    end
    ZeroCrossingRate(i)=sum(TempSignDifference(i,:));
end

% avgZeroCrossingRate Calculation
% Gia ton ypologismo tou avgZeroCrossingRate xrisimopoieitai 
% orthogonio parathiro me mikos oso kai to mikos tou ZeroCrossingRate
% diladi osos einai o arithmos ton frames.
rWin = rectwin(TotalFrameNum);
avgZeroCrossingRate = (1/TotalFrameNum)* (ZeroCrossingRate * rWin);
