function logml = evaluate(data, partition, varargin)
% Calculates the marginal likelihood for the data using a given partition. 
% Author: Paul Blomstedt

if nargin>2
    priorType = 2;
    a = varargin{1};            % set concentration parameter for CRP prior
    if ~(a > 0)
        error('alpha must be > 0')
    end
else
    priorType = 1;                  % the default prior type is uniform over partitions
    a = [];
end

global COUNTS;
global SAMPLEVAR;
global SAMPLEMEAN;
global PARTITION; 
global LOGML_TABLE;
clearGlobalVars;

npops = length(unique(partition));

% Initialize global variables
counts = initialCounts3(partition, data, npops);       
ndat = normalizeData(data); % normalize data for each feature 
[samplevar, samplemean] = initialStats(partition, ndat, npops); 
COUNTS = counts;
SAMPLEVAR = samplevar; 
SAMPLEMEAN = samplemean;  
PARTITION = shiftdim(partition);  
LOGML_TABLE = zeros(npops,1);
clear counts; 
clear partition; 
clear samplevar; 
clear samplemean; 

% compute log-marginal likelihood
updateLogmlTable2(1:npops);
logml = computeTotalLogml2(data,priorType,a);

