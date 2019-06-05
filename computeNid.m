function computeNid(algorithm_type,nTop)
% Evaluate similarity between obtained clusterings using Normalized
% Information Distance.
% Author(s): Paul Blomstedt

load(['Data/MAT/' algorithm_type '_clust_top_' num2str(nTop) '.mat'],'cluster_label');

D_nid = zeros(length(cluster_label));
disp(' ')
for ind1=1:length(cluster_label)
    disp(['Computing NID for experiment no. ' num2str(ind1)])
    for ind2=1:length(cluster_label)
        [~,nid_tmp] = AMI(cluster_label{ind1},cluster_label{ind2});
        D_nid(ind1,ind2) = nid_tmp;
    end
end

save(['Data/MAT/' algorithm_type '_NID_top_' num2str(nTop) '.mat'],'D_nid');
fprintf(['\nSaved results in file Data/MAT/' algorithm_type '_NID_top_' num2str(nTop) '.mat\n\n'])
