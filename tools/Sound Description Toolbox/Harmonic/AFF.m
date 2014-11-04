function f0 = AudioFundamentalFrequencyType(s,standvar,num_frames)

% f0 = AudioFundamentalFrequencyType(s,standvar,num_frames)
% Estimate fundamental frequency
% s is the audiosignal
% standvar contains the parameters of the signal and analysis

% Written By Melanie Jackson
% Version 1.0 5 Feb 2001
% Modified 9 Feb 2001 - Shortened maximum lag and analysis interval size.
% Modified 19th March 2001 - Compatible to variable initialisation

% Example:
% standvar = h_mpeg7init(B,[],[],[],[]);
% f0 = AFF(A,standvar,485);

z = standvar.hopsize;
%%%%%%%%
z = 661;
%%%%%%%%
n = floor(standvar.windowsize);
fs = standvar.fs;
f = [];
A = [];
f0 = [];

maxperiod = 30e-3; % 30 ms
lofreq = 62.5;
hifreq = 1500;

Km = ceil(fs/lofreq); % maximum lag
Kl = floor(fs/hifreq); %minimum lag
% Coping with Startup missing history
K = z; % maxlag for second frame
ProbFrame = ceil(Km/z);

%h = waitbar(0,'Please Wait');
f0 = zeros(num_frames,1);
for frame = 2:num_frames-ProbFrame+1
   m = frame*z;
   den1 = sum(s(m+1:m+n).^2);
   phi = zeros(1,K);
   den =sum(s(m-Kl+2:m-Kl+n).^2);
    for k = Kl:K
        % Normalized Cross Correlation
        den = den-s(m-k+n)^2+s(m-k)^2;
        num =sum(s(m+1:m+n).*s(m-k+1:m-k+n));
        phi(k) = num/(sqrt(den1*den)+eps);
    end
    [mag,index] = max(phi);
    [a,index] = max(phi>0.97*mag);
  
    %index = fs*(index<Kl)+index;

    f0(frame) = fs/index;
    if frame < ProbFrame
        K = K+z;
    elseif frame == ProbFrame
        K = Km;
    end
    %waitbar(frame/num_frames)
end
phi = [];
f0(1) = f0(2); % Use second frame fundamental as 1st frame estimate 
% For last ProbFrame-1 frames use last estimate of fundamental
m = length(s)-n;
den1 = sum(s(m+1:m+n).^2);
for k = Kl:K
    den =sum(s(m-k+1:m-k+n).^2);
    num =sum(s(m+1:m+n).*s(m-k+1:m-k+n));
    phi(k) = num/(sqrt(den1*den)+eps);
end
[mag,index] = max(phi);
%index = (index<Kl)+index;
f0(frame+1:num_frames) = fs/index;
%close(h);
