function meanValues = Mean_SeriesOfScalar(auData, scalingRatio, elementNum, weight_flag, weight, write_flag)

%%% File: Mean_SeriesOfScalar.m
%%
%% --------------  ---------------
%%
%% The function of this subroutine is Series of Mean 
%% of groups of samples. Size of the vector should equal 
%% elementNum, or 0 if Raw is present.
%% 
%% Definition:
%%
%% N = P*Q (?)
%% P -> scaleRatio, Q -> rescaleFactor 
%% MeanValue = \mean_{i=1+(k-1)*N}^{kN} x_i
%% x -> audata
%% If Weight present, ignore samples with zero weight. 
%% If all have zero weight, set to zero by convention. 
%% (N3704)
%%
%%
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

global fid;

%%if nargin < 3, error('constr requires three input arguments'); end
%%if nargin < 5, weight_flag = 0; weight= []; end


%% Initialization
totalSampleNum = length(auData);
len = length(scalingRatio);

meanValues =  [];

if weight_flag == 0
   if len == 1 & scalingRatio(len) == 1& write_flag == 1
      
   else
      lastlenData = 0;
      k = 1;
      len = length(scalingRatio);
      while k <= len
         lenData = lastlenData + scalingRatio(k)*elementNum(k);
         if totalSampleNum < lenData
            N = floor((totalSampleNum - lastlenData)/scalingRatio(k));
            if N >= 0
       			for i = 1:N
         			sn_segment = auData(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                  meanValues = [meanValues mean(sn_segment)];
            	end
            	N1 = lastlenData + N*scalingRatio(k);
            	sn_segment = auData(N1+1:totalSampleNum);
               meanValues = [meanValues, mean(sn_segment)];
            end
         else
            for i = 1:elementNum(k)
              	sn_segment = auData(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
 				   meanValues = [meanValues mean(sn_segment)];
             end
             lastlenData = lenData;
         end

          k=k+1;
      end

	end
else
  	if len == 1 & scalingRatio(len) == 1 & write_flag == 1

   else
      if sum(weight) == 0
         meanValues = [meanValues 0];
   	else
         lastlenData = 0;
         k = 1;
         len = length(scalingRatio);
         while k <= len
          	lenData = lastlenData + scalingRatio(k)*elementNum(k);
         	if totalSampleNum < lenData
            	N = floor((totalSampleNum - lastlenData)/scalingRatio(k));
            	if N >= 0
       				for i = 1:N
                   sn_segment = auData(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                   weight_segment = weight(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                   sn_weight = weight_segment.*sn_segment;
            	    mean_sn_weight = sum(sn_weight)/sum(weight_segment);
						 meanValues = [meanValues mean_sn_weight];
            		end
               	N1 = lastlenData + N*scalingRatio(k);
                  sn_segment = auData(N1+1:totalSampleNum);
                  weight_segment = weight(N1+1:totalSampleNum);
                  sn_weight = weight_segment.*sn_segment;
                  mean_sn_weight = sum(sn_weight)/sum(weight_segment);
						meanValues = [meanValues mean_sn_weight];
            	end
         	else
            	for i = 1:elementNum(k)
                  sn_segment = auData(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                  weight_segment =weight(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                  mean_sn_weight = sum(sn_weight)/sum(weight_segment);
                  meanValues = [meanValues mean_sn_weight];
               end
               lastlenData = lenData;
            end

          	k=k+1;
       	end
       end;

    end
 end
 
 
    
       

   
