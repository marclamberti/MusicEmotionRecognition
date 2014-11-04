function y = ma_test_create_wav(fs,freq,loudness,duration,filename)
%%
%% y = ma_test_create_wav(fs,freq,loudness,duration,filename)
%%
%% create a wav file for testing some of the functions
%% sinus with different frequencies and loudness
%% default settings: frequency based on bark scale
%%
%% usage e.g.:
%%            sound(ma_test_create_wav) %% listen to defaults
%%            freq = repmat([150,350,350,0],1,10);
%%            sound(ma_test_create_wav(22050,freq,.2,0.2)) %% periodic pattern
%%
%% fs       ... sampling frequency (default = 22050 Hz)
%% freq     ... frequencies of sinus tones [Hz] (use [] for default)
%% loudness ... loudness levels to use [0..1] (use [] for default)
%% duration ... duration of sinus tones in sec (use [] for default = 0.2 sec)
%% filename ... if given, wav file is created
%%
%% y        ... output vector [0..1] (e.g. "sound(y,fs)")

%% elias 16.5.2004

if nargout==0 & nargin==0,
    disp('testing: ma_test_create_testwav (listen to defaults)')
    sound(ma_test_create_wav);
    return
end

if nargin<2 | isempty(freq),
    freq = [100 200 300 400 510 630 770 920 1080 1270 1480 ...
            1720 2000 2320 2700 3150 3700 4400 5300 6400 ...
            7700 9500 12000 15500]; %% see zwicker & fastl (1999)
end

if nargin<3 | isempty(loudness),
    loudness = linspace(0.2,1,3);
end
    
if nargin<4 | isempty(duration),
    duration = .2;
end

if nargin==0,
    fs = 22050;
end

freq = freq(freq<fs/2);

y = zeros(length(freq)*length(loudness)*round(duration*fs),1);

w = 0.5 * (1 - cos(2*pi*(0:1024-1)'/(1024-1))); %% hann

idx = 1:round(fs*duration);
for i=1:length(loudness),
    for j=1:length(freq),
        y_tmp = sin(linspace(0,freq(j)*2*pi*duration,round(fs*duration)));
        y_tmp(1:512) = y_tmp(1:512).*w(1:512)';
        y_tmp(end-512+1:end) = y_tmp(end-512+1:end).*w(513:end)';
        y(idx) = y_tmp*loudness(i);
        idx = idx + round(fs*duration);
    end
end

if nargin==5,
    wavwrite(y,fs,16,filename)
end