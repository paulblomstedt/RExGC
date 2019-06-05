function createData

% This function downloads the data, creates a folder structure and selects the 
% A-AFFY-44 experiments to be used in the analyses and creates ground truth 
% matrices based on EFOs. The selection of experiments is done in two stages: 
%
% In the first stage all experiments for which both measurement data and experiment 
% analytics are identified. Further, experiments with a very small number of genes 
% are discarded (most experiments have > 54670 genes, so this is set as the limit).
%
% In the second stage the selection is based on availability of EFOs to be used as 
% ground truth. Here, the selection criteria are 
% 1. presence of EFO values corresponding to at least one of: cell type, disease, organism part
% 2. unique value for a given EFO type
% 3. each EFO value present in at least two experiments
%
% The output for the three EFO types are stored in a separate files: 
% 'Data/MAT/G44_cell_type.mat', 'Data/MAT/G44_disease.mat', 'Data/MAT/G44_organism_part.mat', each of 
% which contains a binary similarity matrix to be used as ground truth, as well as a list of accession 
% numbers of the experiments which satisfy the selection criteria. A combined ground truth matrix is 
% also generated and stored.
%
% Author(s): Paul Blomstedt, Ritabrata Dutta

%% Download data and create folder structure

system('Rscript downloadData.r')


%% Initial selection of experiments 

fprintf('\nInitial selection of experiments in process, this may take a while...\n\n')

% Load data
tabFilesRaw=dir('Data/TAB/Raw/*.tab');                        % files of gene expression measurement data
tabFilesAnalysed=dir('Data/TAB/Analysed/*.tab');            % files of experiment analytics (t-statistics and p-values)

% Extract accession numbers from files of measurement data
expRaw = cell(length(tabFilesRaw),1);
for k = 1:length(tabFilesRaw)
    a = strsplit(tabFilesRaw(k).name,'-');
    expRaw(k) = {[char(a(2)) '-' char(a(3)) '-' char(a(4))]};
end

% Extract accession numbers from files of experiment analytics
expAnalysed = cell(length(tabFilesAnalysed),1);
for k = 1:length(tabFilesAnalysed)
    a = strsplit(tabFilesAnalysed(k).name,'-');
    expAnalysed(k) = {[char(a(2)) '-' char(a(3)) '-' char(a(4))]};
end

% Select experiments if both measurement data and experiment analytics exist
[~,indRaw] = intersect(expRaw,expAnalysed); 

% Select experiments with > 54670 genes
nGenes = zeros(1,length(tabFilesRaw));                                         % no. of genes contained in each experiment
for k=1:length(tabFilesRaw)
    k
    myData = importdata(['Data/TAB/Raw/' tabFilesRaw(k).name]);
    nGenes(k) = length(myData.textdata(2:end,1));  
end
indGenes = find(nGenes > 54670);

% List accession numbers of experiments which satisfy the above criteria
indSelect = intersect(indRaw,indGenes);  
expNamesGE = expRaw(indSelect);         

%% Selection of experiments and creation of ground truth based on EFOs

fprintf('\nConstructing the ground truth...\n\n')

EFOtypes_groundtr = {'cell_type','disease','organism_part'};
G_combined = zeros(length(expNamesGE));
indCombined = [];
for ind = 1:length(EFOtypes_groundtr)
    
    % Load EFO data 
    EFOtype = EFOtypes_groundtr{ind};
    efoFile = readtable(['Data/EFO/' EFOtype '.csv']);
    
    % 1. Select experiments for which the given EFO type is present
    [expNamesEFO,indEFO,n] = unique(efoFile(:,1)); 
    annotation = efoFile(indEFO,2);
    [~,indEFO_exists] = intersect(table2cell(expNamesEFO),expNamesGE);
    
    % 2. Select experiments with a unique value for a given EFO type 
    nEFO = accumarray(n, 1, [], @sum); 
    indEFO_unique = find(nEFO==1);  
    
    % Find subset of experiments by combining steps 1 and 2.
    indSubset = intersect(indEFO_exists,indEFO_unique); 
    expNamesEFO = table2cell(expNamesEFO(indSubset,:));
    annotation = table2cell(annotation(indSubset,:));
    
    % 3. Select experiments with the same EFO value present in at least two experiments 
    cmpEFO = zeros(length(annotation));     
    for i = 1:length(annotation)
        for j = 1:length(annotation)
            cmpEFO(i,j) = strcmp(annotation{i},annotation{j});
        end
    end
    indEFO_repeat = find(sum(cmpEFO)>1);
    expNamesEFO = expNamesEFO(indEFO_repeat);
    annotation = annotation(indEFO_repeat);
    
    % Creation of ground truth
    [~,indGE] = intersect(expNamesGE,expNamesEFO);
    G = zeros(length(annotation));     
    for i = 1:length(annotation)
        for j = 1:length(annotation)
            G(i,j) = strcmp(annotation{i},annotation{j});
        end
    end

    G_combined(indGE,indGE) = G_combined(indGE,indGE)+G;    % combined similarity matrix for all EFO types
    indCombined = union(indCombined,indGE);                 % combined index set for final selection of experiments
    
    % Assign results to variables
    eval(['G_' EFOtype ' =  G;']);
    eval(['indGE_' EFOtype ' =  indGE;']);
    eval(['annotation_' EFOtype ' =  annotation;']);
    eval(['expNames_' EFOtype ' =  expNamesEFO;']);
end


%% Discard experiments which are not used and save 

G = G_combined(indCombined,indCombined);          % final combined similarity matrix
expNames = expNamesGE(indCombined);               % final list of experiments
save('Data/MAT/G44_EFO_combined.mat','G','expNames');
fprintf('\nSaved combined ground truth matrix in Data/MAT/G44_EFO_combined.mat\n\n')

for ind = 1:length(EFOtypes_groundtr)
    EFOtype = EFOtypes_groundtr{ind};
    eval(['G =  G_' EFOtype ';']);
    eval(['expNamesEFO = expNames_' EFOtype ';']);
    eval(['annotation = annotation_' EFOtype ';']);
    [~,indKeep] = intersect(expNames,expNamesEFO);
    save(['Data/MAT/G44_' EFOtype '.mat'],'G','expNames','annotation','indKeep');
    fprintf(['\nSaved EFO-specific ground truth matrix in  Data/MAT/G44_' EFOtype '.mat\n\n'])
end


