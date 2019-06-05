
%% 1. Cluster the data

% rexClustering('PPM',5);               % precomputed output available    
rexClustering('k-means',5);

%% 2a. Compute similarity matrices for clusterings using Normalized Information Distance

% computeNid('PPM',5);                  % precomputed output available
computeNid('k-means',5);

%% 2b. Compute similarity matrices for baseline approaches 

% evaluateLogML;                        % precomputed output available
% DEpValCorrelation;                    % requires analytics data to be downloaded (see step 1 in README.txt), precomputed output available

%% 3. Visualize results

%  algorithm = {'PPM','LogML','DE'};    % reproduces Fig. 1 of the paper, the required input is available precomputed 
algorithm = {'k-means','LogML','DE'};

EFO = {'cell_type', 'disease', 'organism_part'};
legStr = [{'Random'},algorithm];

for i = 1:length(EFO)
    rexVisualization(EFO{i},algorithm,5)
    legend(legStr,'Location','best')
end