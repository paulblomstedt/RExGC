function evaluateLogML
% Similarity matrix for likelihood-based approach
% Author(s): Paul Blomstedt

load('Data/MAT/PPM_clust_top_5.mat');
load('Data/MAT/44_top_5.mat','dataProcessed');

nExp = length(cluster_label);
D_lml = zeros(nExp);
for i = 1:nExp
    disp(['Evaluating log-likelihood for experiment no. ' num2str(i)])
    for j = 1:nExp
        D_lml(i,j) = evaluate(dataProcessed{i},cluster_label{j});
    end
end

maxVal= max(max(D_lml));
D_lml(logical(eye(nExp))) = maxVal+1; % make sure the diagonal has the max value

save('Data/MAT/LogML.mat','D_lml');
fprintf('\nSaved results in file Data/MAT/LogML.mat\n\n')
