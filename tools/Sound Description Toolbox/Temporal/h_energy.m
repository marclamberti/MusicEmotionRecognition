% function [energy_bp] = Fcalculenv(data_v, sr_hz, cutfreq_hz, dsfact)
%
% INPUTS
% ======
% - data_v      : vector containing the data of the soundfile
% - sr_hz       : sampling rate of the soundfile
% - cutfreq_hz  : cutting frequency for low-pass filtering of the energy
% - dsfact      : down-sampling factor for the energy [integer] (1=Fe, 2=Fe/2, 3=Fe/3, ...)
% 
% OUTPUTS
% =======
% - energy_bp   : breakpoint function 
%                 [first colum: time [second] | second column: energy value]
%
% Target:   MP7-XM version
% Author:   CUIDADO/IRCAM/ G. Peeters 
% LastEdit: 2001/03/12
%

function [energy_bp] = Fcalculenv(data_v, sr_hz, cutfreq_hz, dsfact)

  L_n    = round(sr_hz / (2*cutfreq_hz));
  L_n    = L_n + ~rem(L_n,2);
  LD_n   = (L_n-1)/2;
  STEP_n = dsfact;
  
  mark_n_v   = [1+LD_n : STEP_n : length(data_v)-LD_n];
  time_sec_v = (mark_n_v-1)/sr_hz;
  
  nb_frames = length(mark_n_v);
  energy_v  = zeros(1, nb_frames);
  	
  for index = 1:nb_frames
    
    n = mark_n_v(index);
    
    signal    = data_v(n-LD_n:n+LD_n);
    signal_dc = signal - mean(signal);
    energy_v(index) = sqrt( sum(signal_dc.^2) / L_n );
      
  end
  
  energy_bp = [time_sec_v(:), energy_v(:)];