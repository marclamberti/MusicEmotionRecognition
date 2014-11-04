% Sinartisi save_feature_matrix: Kalei tin sinartisi feature_matrix gia
% na swsei ton pinaka features enos arxeiou ixou ston disko se ena arxeio
% *.sfm.

function [] = save_feature_matrix(filename)

fid = fopen([filename '.fm'],'w+');
Z = ComputeFeatureMatrix(filename);
fprintf(fid, '%4.8f\n',Z);
fclose(fid);

clear Z;