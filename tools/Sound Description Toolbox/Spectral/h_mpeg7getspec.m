% Modified 30/04/2002 by Thorsten Kastner - zero pad (length overlap) at start and end for overlap removed; 
%                                           because: there's no overlap in ASF
function  [fftout,phase] = mpeg7getspec(data,v)
% data = data(1:5004);
if size(data,1)==1
    data = data';
end

hs = mean(v.hopsize);
hops = length(v.hopsize);
num_f = ceil(length(data)/hs); 
rem_f = mod(num_f,hops);
pad = sum(v.hopsize(1:rem_f)) - (length(data)-(num_f-rem_f)*hs);
data = [data; zeros(pad,1)];
num_samples = length(data);
fstart = 1;
fftout = [];
NormWindow = sum(v.window.*v.window);
spec = (h_specgram2(data,v.FFTsize,v.fs,v.window,v.hopsize))/sqrt(v.FFTsize*NormWindow);
fftout = abs(spec);
phase = angle(spec);
%  figure(30); imagesc([1 size(fftout(2))],[0 v.fs/2],fftout)
% Need to compensate first and last frames for the zero padding

if pad
    fftout(:,end) = fftout(:,end)*sqrt(v.windowsize/(v.windowsize-pad));
end
