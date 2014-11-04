function featureMatrix = ComputeFeatureMatrix(filename)

% Initialize
startup;
format('long');
warning('off');
[A,B,C] = wavread(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AP
AudioPower_SeriesOfScalar = AP(A,size(A),B,4,64,[]);
meanAP = mean(mean(AudioPower_SeriesOfScalar,2));
covAP = cov(mean(AudioPower_SeriesOfScalar,2));
meanDiffAP = mean(diff(mean(AudioPower_SeriesOfScalar,2)));
covDiffAP = cov(diff(mean(AudioPower_SeriesOfScalar,2)));

clear AudioPower_SeriesOfScalar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harmonic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AFF
standvar = h_mpeg7init(B,[],[],[],[]);
f0 = AFF(A,standvar,485);
meanf0 = mean(f0);
covf0 = cov(f0);
meanDifff0 = mean(diff(f0));
covDifff0 = cov(diff(f0));

clear f0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perceptual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TL_SONE
p = struct('fs',B,'do_visu',0,'do_sone',1,'do_spread',1,'outerear','terhardt','fft_size',256,'hopsize',128,'bark_type','table','dB_max',96);
[sone, Ntot, p] = TL_SONE(A,p);
meanTL = mean(Ntot);
covTL = cov(Ntot);
meanDiffTL = mean(diff(Ntot));
covDiffTL = cov(diff(Ntot));
meanSONE = mean(sone,2); meanSONE = meanSONE(1:8);
for i=1:8 covSONE(i) = cov(sone(i,:)); end
for i=1:8 meanDiffSONE(i) = mean(diff(sone(i,:))); end
for i=1:8 covDiffSONE(i) = cov(diff(sone(i,:))); end

clear Ntot;
clear sone;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ASC
AudioSpectrumCentroid = ASC(filename,'PT10N1000F',0,[]);
meanASC = mean(AudioSpectrumCentroid);
covASC = cov(AudioSpectrumCentroid);
meanDiffASC = mean(diff(AudioSpectrumCentroid));
covDiffASC = cov(diff(AudioSpectrumCentroid));

clear AudioSpectrumCentroid;

% ASR
AudioSpectrumRolloff = ASR(filename,0.016);
meanASR = mean(AudioSpectrumRolloff);
covASR = cov(AudioSpectrumRolloff);
meanDiffASR = mean(diff(AudioSpectrumRolloff));
covDiffASR = cov(diff(AudioSpectrumRolloff));

clear AudioSpectrumRolloff;

% ASS
AudioSpectrumSpread = ASS(filename,'PT10N1000F',0,[]);
meanASS = mean(AudioSpectrumSpread);
covASS = cov(AudioSpectrumSpread);
meanDiffASS = mean(diff(AudioSpectrumSpread));
covDiffASS = cov(diff(AudioSpectrumSpread));

clear AudioSpectrumSpread;

% MFCC
[ceps,freqresp,fb,fbrecon,freqrecon] = MFCC(A, B);
[C1 C2] = size(ceps);
for i=1:C1 for j=1:C2 if(isinf(ceps(i,j)) || isnan(ceps(i,j))) ceps(i,j) = 0; end; end; end;
meanMFCC = mean(ceps,2);
for i=1:C1 covMFCC(i) = cov(ceps(i,:)); end
for i=1:C1 meanDiffMFCC(i) = mean(diff(ceps(i,:))); end
for i=1:C1 covDiffMFCC(i) = cov(diff(ceps(i,:))); end

clear ceps;
clear freqresp;
clear fb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AC
ACcoefs = AC(A);

% LAT
energy_bp = h_energy(A, B, 2, 2000);
LogAttackTime = LAT(energy_bp, 50);

% TC
TemporalCentroid = TC(energy_bp);

clear energy_bp;

% ZCR
[ZCR1 avZCR] = ZCR(A,0.016);
meanZCR = mean(ZCR1);
covZCR = cov(ZCR1);
meanDiffZCR = mean(diff(ZCR1));
covDiffZCR = cov(diff(ZCR1));

clear ZCR1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Various %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ASF
[AudioSpectrumFlatness ,lo_edge, hi_edge, XMLFile ] = ASF(A,B,'PT10N1000F',250,500,0,[]);
[ASF1 ASF2] = size(AudioSpectrumFlatness);
meanASF = mean(AudioSpectrumFlatness,1);
for i=1:ASF2 covASF(i) = cov(AudioSpectrumFlatness(i,:)); end
for i=1:ASF2 meanDiffASF(i) = mean(diff(AudioSpectrumFlatness(i,:))); end
for i=1:ASF2 covDiffASF(i) = cov(diff(AudioSpectrumFlatness(i,:))); end

clear AudioSpectrumFlatness;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put all features together %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

featureMatrix = [];
featureMatrix = [featureMatrix meanAP covAP meanDiffAP covDiffAP];
featureMatrix = [featureMatrix meanf0 covf0 meanDifff0 covDifff0];
featureMatrix = [featureMatrix meanTL covTL meanDiffTL covDiffTL];
featureMatrix = [featureMatrix meanSONE' covSONE meanDiffSONE covDiffSONE];
featureMatrix = [featureMatrix meanASC covASC meanDiffASC covDiffASC];
featureMatrix = [featureMatrix meanASR covASR meanDiffASR covDiffASR];
featureMatrix = [featureMatrix meanASS covASS meanDiffASS covDiffASS];
featureMatrix = [featureMatrix meanMFCC' covMFCC meanDiffMFCC covDiffMFCC];
featureMatrix = [featureMatrix ACcoefs'];
featureMatrix = [featureMatrix LogAttackTime];
featureMatrix = [featureMatrix TemporalCentroid];
featureMatrix = [featureMatrix meanZCR covZCR meanDiffZCR covDiffZCR];
featureMatrix = [featureMatrix meanASF covASF meanDiffASF covDiffASF];

clear meanAP covAP meanDiffAP covDiffAP;
clear meanf0 covf0 meanDifff0 covDifff0;
clear meanTL covTL meanDiffTL covDiffTL;
clear meanSONE covSONE meanDiffSONE covDiffSONE;
clear meanASC covASC meanDiffASC covDiffASC;
clear meanASR covASR meanDiffASR covDiffASR;
clear meanASS covASS meanDiffASS covDiffASS;
clear meanMFCC covMFCC meanDiffMFCC covDiffMFCC;
clear ACcoefs LogAttackTime TemporalCentroid;
clear meanZCR covZCR meanDiffZCR covDiffZCR;
clear meanASF covASF meanDiffASF covDiffASF;
