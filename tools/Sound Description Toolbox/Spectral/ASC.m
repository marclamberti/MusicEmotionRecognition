function [AudioSpectrumCentroid, XMLFile] =AudioSpectrumCentroidType(audioFile,hopSize,writeXML,XMLFile)


% Example: 
% AudioSpectrumCentroid = ASC('BaCl.mf.C4B4_3.wav','PT10N1000F',0,[]);

% This function describes the centre of gravity of log-frequency power spectrum
% Audiosignal is Audio data.

% Written by Melanie Jackson
% Version 1.0 12 Jan 2001
% Modified 16 March 2001 - Removed spectrum extraction to generic function
% Modified 16/04/2002 by Thibaut Sacreste - add XML generation
% Modified 19/04/2002 by Thibaut Sacreste - changed to be a stand alone function
% Modified 03/05/2003 by Holger Crysandt - power spectrum bug-fix, see w5048 for details

%--------------------------------------------------------------------
% audioFile is the name of the audio file to process
% 2 types of files can be read: .wav and .au 
% writeXML is a flag for the generation of the XML file
% writeXML=0 -> no generation
% writeXML=1 -> generation
% XMLFile is the name of the XML file to be generated (optional)
%--------------------------------------------------------------------


% Initialisation:
if nargin<4 writeXML=0;end
if nargin<2 hopSize='PT10N1000F'; end;

% Read in audio file
if audioFile(end-3:end)=='.wav'
    [audioData, sr] = wavread(audioFile);
elseif audioFile(end-2:end)=='.au'
    [audioData, sr] = auread(audioFile);
else
    Error = 'Incorrect filename'
    return
end

% Descriptors only deal with monaural recordings 
if size(audioData,2)>1 
    % in the wavread function the second dimension contains the number of channels
    audioData = mean(audioData')';
end

% Start calculating descriptors

% Hopsize conversion
% format PT10N1000F

i = find(hopSize=='N');
hop = [str2num(hopSize(3:i-1)) str2num(hopSize(i+1:end-1))];


standvar = h_mpeg7init(sr,hop);


%------------------------------------------------------------------
% STFT with 1/3 overlap and window size of three times the hopsize.
% Zero padding of the last few frames will occur, to ensure there is one spectral frame
% for each corresponding power estimate in the power descriptor. Zero padding will also occur
% at the start for the same purpose, since 1/3 overlap is used, the emphasis of the
% information lay in the centre of the window.

[fftout,phase] = h_mpeg7getspec(audioData,standvar);


%------------------------------------------------------------------
% AudioSpectrumCentroid calculation
loedge = 62.5;
numframes = size(fftout,2); 
N = standvar.FFTsize;
fs = standvar.fs;
fft_freq =(0:N/2)*fs/N;
% replace frequencies less than 62.5Hz with nominal freq 31.25Hz
num_less_loedge = sum(fft_freq<loedge);
fft_freq = [31.25 fft_freq(num_less_loedge+1:end)];
% determine log scaled frequencies relative to 1kHz
fft_freq_log = log2(fft_freq/1000);

powers = fftout.^2;
powers(2:end-1,:) = 2*powers(2:end-1,:);
if num_less_loedge > 1
        powers = [sum(powers(1:num_less_loedge,:)); powers(num_less_loedge+1:end,:)];
end

AudioSpectrumCentroid = sum((fft_freq_log'*ones(1,numframes)).*powers)./(sum(powers)+eps);


%---------------------
%XML generation:

if writeXML
    if ~exist('XMLFile')
        XMLFile=h_ASCtoXML(AudioSpectrumCentroid,hopSize);
    else
        XMLFile=h_ASCtoXML(AudioSpectrumCentroid,hopSize,XMLFile);
    end 
end    

