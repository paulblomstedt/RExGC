function rexClustering(algorithm_type,nTop)
% Clustering of gene expression experiments. 
% Author(s): Paul Blomstedt, Ritabrata Dutta 

% Load data and initialize output variables
load(['Data/MAT/44_top_' num2str(nTop) '.mat'],'dataProcessed');

Y = dataProcessed; clear('dataProcessed');
cluster_label = cell(1,length(Y));


% Perform clustering on each experiment separately
nclust = ceil(sqrt(size(Y{1},1))/2);            % for k-means, set fixed number of clusters heuristically 
range = 2:ceil(sqrt(size(Y{1},1)));             % for PPM on reduced model space, set range for number of clusters  
solutions = zeros(size(Y{1},1),length(range));
logML = zeros(1,length(range));

t = 0;                                  % set counter
disp(' ')
for index = 1:length(Y)
    t = t+1;
    disp(['Clustering experiment no. ' num2str(t) ' using ' algorithm_type])
    switch algorithm_type
        case 'k-means'
            %warning off stats:kmeans:FailedToConverge 
            cluster_label{1,index} = kmeans(Y{1,index},nclust,'Distance','sqeuclidean');
        case 'PPMkm'
            %warning off stats:kmeans:FailedToConverge 
            for i=1:length(range)
                solutions(:,i) = kmeans(Y{1,index},range(i),'Distance','sqeuclidean');
                logML(i) = evaluate(Y{1,index},solutions(:,i),1);
            end
            [~,I] = max(logML);
            cluster_label{1,index} = solutions(:,I);
        case 'PPM'
            cluster_label{1,index} = analyzeData(Y{1,index},range(end));
    end
end

%% Save output

save(['Data/MAT/' algorithm_type '_clust_top_' num2str(nTop) '.mat'],'cluster_label');
fprintf (['\nSaved results in file Data/MAT/' algorithm_type '_clust_top_' num2str(nTop) '.mat\n\n'])



