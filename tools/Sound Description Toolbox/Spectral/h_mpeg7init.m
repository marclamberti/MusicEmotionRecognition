function standvar = mpeg7init(fs,hopsize,windowsize,window,FFTsize)

% This function creates a structure of the default values
% to be used throughout the descriptors

% Written by Melanie Jackson

% Version 1 15th March 2001
% Modified 30/04/2002 by Thorsten Kastner - changed standard hopsize to 30ms
%                                         - set hopsize = windowsize (no overlap for ASF)
  
if nargin == 0
    Error = 'need to specify the sampling frequency: standvar = mpeg7init(fs)'
    standvar = [];
    return 
end

try	% variable defined
    if isempty(hopsize)			% variable defined but no value
      hopsize = [fs*30 1000]; % Changed hopsize: 30ms recommend
    else
        hopsize = [fs*hopsize(1) hopsize(2)];
    end   
catch
    hopsize = [fs*30 1000] ;
end
[q, n, d] = h_fraction(hopsize(1),hopsize(2));
% This next section of code is to ensure minimal stray from the sampling period
% This is done by interleaving minor hopsize with major
% e.g if  10 10 10 10 10 11 11 are the hops then the pattern should be
% 10 10 10 11 10 10 11
if n==0
    hopsize = q;
elseif d-n>=n
    k = ceil((d-n)/n); % ratio of q to q+1 occurence
    fr = floor((d-n)/(k)); % number of sub sequences required
    frr = [q*ones(1,k) q+1]; % subsequence - length k+1
    hopsize = reshape(frr'*ones(1,fr),1,(k+1)*fr); % append subsequences to each other
    hopsize = [hopsize q*ones(1,d-n-k*fr) (q+1)*ones(1,n-fr)]; % attach extra values
else  
    k = ceil(n/(d-n));
    fr = floor(n/k);
    frr = [(q+1)*ones(1,k) q];
    hopsize = reshape(frr'*ones(1,fr),1,(k+1)*fr);
    hopsize = [hopsize (q+1)*ones(1,n-k*fr) q*ones(1,d-n-fr)];
end


try
  if isempty(windowsize) 	% variable defined but no value
      windowsize = ceil(fs*30/1000); % Changed:Windowsize fixed at 30ms
    end
catch
  windowsize = ceil(fs*30/1000); % Changed: Windowsize fixed at 30ms
end

try
  if isempty(FFTsize)			% variable defined but no value
    FFTsize = 2^nextpow2(windowsize);
  end
catch
  FFTsize = 2^nextpow2(windowsize);
end

try
  if isempty(window)			% variable defined but no value
    window = [hamming(windowsize)];
  end
catch
  window = [hamming(windowsize)];
end

standvar = struct('fs',fs,'hopsize',hopsize,'windowsize',windowsize,'window',window,'FFTsize',FFTsize);
return
