function [AudioSpectrumFlatness ,lo_edge, hi_edge, XMLFile ] = AudioSpectrumFlatnessType(audioData, sr,hopSize,loEdge,hiEdge,writeXML,XMLFile)


% Example:
% [AudioSpectrumFlatness ,lo_edge, hi_edge, XMLFile ] = ASF('BaCl.mf.C4B4_3.wav','PT10N1000F',250,500,0,[]);
% Recommended low and high-edge: (250,500) (500,1000) (1000,2000)
% (2000,4000)

% For tonal signals: ASF -> 0
% For noisy signal: ASF -> 1

% This function describes the spectral flatness measure of the audio signal
% The frequency range is divided into numbands logarithmically spaced bands
% of 1/4 octave width, starting at loEdge up to hiEdge

% Written by Melanie Jackson & Juergen Herre
% Version 2.0 25 July 2001

% Modified by Melanie Jackson 10 December 2001
% Modified by MJ 17 January 2002
% Modified 16/04/2002 by Thibaut Sacreste - add XML generation
% Modified 19/04/2002 by Thibaut Sacreste - changed to be a stand alone function
% Modified 30/04/2002 by Thorsten Kastner - changed function call (loEdge and hiEdge can be set)
%                                         - modified check for loEdge and hiEdge
% Modified 12/06/2002 by Thorsten Kastner - adapted coefficient grouping to variable loEdge and hiEdge
%                                         - values for loEdge and hiEdge will be recalculated if necessary

%--------------------------------------------------------------------
% audioFile is the name of the audio file to process
% 2 types of files can be read: .wav and .au 
% writeXML is a flag for the generation of the XML file
% writeXML=0 -> no generation
% writeXML=1 -> generation
% XMLFile is the name of the XML file to be generated (optional)
%--------------------------------------------------------------------
% Initialisation:


% Read in audio file
% if audioFile(end-3:end)=='.wav'
%     [audioData, sr] = wavread(audioFile);
% elseif audioFile(end-2:end)=='.au'
%     [audioData, sr] = auread(audioFile);
% else
%     Error = 'Incorrect filename'
%     return
% end

% Descriptors only deal with monaural recordings 
if size(audioData,2)>1 
    % in the wavread function the second dimension contains the number of channels
    audioData = mean(audioData')';
end

% Start calculating descriptors

% Hopsize conversion
% format PT10N1000F
hop = '';
if (hopSize)
  i = find(hopSize=='N');
  hop = [str2num(hopSize(3:i-1)) str2num(hopSize(i+1:end-1))];
end


standvar = h_mpeg7init(sr,hop);

%Check loEdge and hiEdge

if isempty(loEdge)	% variable defined but no value
  loEdge = 250;
else
  loEdge = 2^(floor(log2(loEdge/1000)*4)/4)*1000 ; % Setting exact value for loEdge, rounding to next lower frequency if necessary
  if (loEdge < 250 ) 
    loEdge = 250;       % No extraction below 250Hz
  end
end 


if isempty(hiEdge)	% Variable defined but no value
  hiEdge = 16000;       % Setting default hiedge
end
if (hiEdge >= sr/2)     % Check if value for hiedge is valid
  hiEdge =  min(hiEdge,sr/2-1); % sr/2-1 : Skipping extraction up to sr/2; Not possible due to 5 percent band overlap
end
hiEdge = 2^(floor(log2(hiEdge/1000)*4)/4)*1000 ; % Setting exact value for hiEdge, rounding to next lower frequency if necessary
if (hiEdge*1.05 >= sr/2) hiEdge = hiEdge/2^0.25 ; end; %Now it's possible to check if hiEdge is valid
						       %If hiedge plus band overlap greater than sr/2 skip highest band

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%------------------------------------------------------------------
% STFT with no overlap; window size equal to  hopsize.
% Zero padding of the last few frames will occur, to ensure there is one spectral frame
% for each corresponding power estimate in the power descriptor. 

[fftout,phase] = h_mpeg7getspec(audioData,standvar);

%-----------------------------------
% AudioSpectrumFlatness calculation:

fs = standvar.fs;
N = standvar.FFTsize;

numbands = floor(4*log2(hiEdge/loEdge));
firstband = round(log2(loEdge/1000)*4);
overlap = 0.05;

grpsize = 1;
for k = 1:numbands
  f_lo = loEdge * (2^((k-1)/4)) * (1-overlap);
  f_hi = loEdge * (2^((k  )/4)) * (1+overlap);
  i_lo = round( f_lo/(fs/N) ) + 1;
  i_hi = round( f_hi/(fs/N) ) + 1;

  % Rounding of upper index according due to coefficient grouping
  if (k+firstband-1 >= 0)                   %Start grouping at 1kHz 
    grpsize = 2^ceil( (k+firstband )/4);
    i_hi = round((i_hi-i_lo+1)/grpsize)*grpsize + i_lo-1 ;
  else
    grpsize = 1;
  end
  tmp = fftout(i_lo:i_hi,:) .^ 2;         % PSD coefficients
  ncoeffs = i_hi - i_lo + 1;
 
  if (k+firstband-1 >= 0)                   % Coefficient grouping
    tmp2 = tmp(1:grpsize:ncoeffs,:);
    for g=2:grpsize
      tmp2 = tmp2 + tmp(g:grpsize:ncoeffs,:) ;
    end
    tmp = tmp2;
  end
  % Actual calculation
  ncoeffs = ncoeffs/grpsize ;
  tmp = tmp + 1e-50;       % avoid underflow for zero signals
  gm(k,:) = exp( sum(log(tmp))/ncoeffs ); % log processing avoids overflow
  am(k,:) = sum(tmp) / ncoeffs;
end
AudioSpectrumFlatness = (gm./am)';
lo_edge = loEdge;
hi_edge = hiEdge;

%---------------------
%XML generation:

if writeXML
    if ~exist('XMLFile')
      
        XMLFile=h_ASFtoXML(AudioSpectrumFlatness',loEdge,hiEdge,hopSize); 
    else
      XMLFile
        XMLFile=h_ASFtoXML(AudioSpectrumFlatness',loEdge,hiEdge,hopSize,XMLFile);
    end
end
