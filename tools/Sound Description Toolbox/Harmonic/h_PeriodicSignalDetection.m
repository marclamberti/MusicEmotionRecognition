function [minj_pos, H_k] = h_PeriodicSignalDetection(signal,samplingRate) 

%% File: h_PeriodicSignalDetection.m
%% 
%% The algorithm is:
%% a)	Calculate:
%%       r_k(j) = sum(s(i+j)-s(i)).^2)/(sum(s(i).^2)+sum(s(i+j).^2))
%% b)	Choose the minimum, and subtract it from 1:
%%        H_k = 1- min(r(j))
%% This value is 1 for a purely periodic signal, and it should be 
%% close to 0 for white noise.
%% The integration window N should be chosen equal to the largest 
%% expected period (by default: 40 ms = 1/25Hz). The estimate can 
%% be refined by replacing each local minimum of r_k(j)  by the minimum 
%% of a 3-point parabolic fit centered upon it.
%%
%%
%% input:
%% -signal:                          incoming signal
%% -samplingRate:                    sampling Rate
%% output:
%% -minj_pos:                        the minimum position for j
%% -H_k:                             H_k == 1 purely periodic signal
%%                                   H_k == 0 white noise
%%                                   H_k ~= 1 & H_k ~=0 harmonic signal
%%
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

%% BUILD THE TIME VECTOR
sampleNum = length(signal);
time = [0:sampleNum-1]*(1/samplingRate);
min_pos = [];
min_rk = [];



%% === Calculate: r(j) = sum(s(i+j)-s(i)).^2)/(sum(s(i).^2)+sum(s(i+j).^2))===
total_power = sum(signal.^2);
for k = 1:sampleNum-1
	part_power = sum(signal(k+1:sampleNum).^2);
    %%%%%%%%%%%%% MY CODE!!!!!!!!!!!
    if(total_power+part_power == 0) 
       r(k) = 0;
       break;
    end 
    
   r(k)= sum((signal(k+1:sampleNum)-signal(1:sampleNum-k)).^2)/(total_power+part_power);
end
for k = 2:sampleNum-20
   if (((r(k)-r(k-1)) < 0) & ((r(k)-r(k+1)) <0))
      %% PARABOLIC APPROXIMATION
   	 a = (((r(k-1)-r(k))/(time(k-1)-time(k))) - ...
         ((r(k-1)-r(k+1))/(time(k-1)-time(k+1))))/(time(k)-time(k+1));
      b = ((r(k-1)-r(k+1))/(time(k-1)-time(k+1))) - a*(time(k-1)-time(k+1));
      min_position = (-b/(2*a)+time(k))*samplingRate;
      estimated_r = interp1(r,min_position,'linear');
   %   min_position = k;
      min_pos = [min_pos min_position];
   %   min_rk  = [min_rk r(k)];
      min_rk = [min_rk estimated_r];
   end
end
[mrk, mpos] = min(min_rk);
minj_pos = min_pos(mpos);
H_k = 1 - mrk;
%if H_k <= 1.e-4
%   sprintf('%s','Signal is white noise')
%elseif H_k == 1
%   sprintf('%s','Signal is a purely periodic signal')
%else
%   sprintf('%s','Signal is a harmonic signal')
%end

