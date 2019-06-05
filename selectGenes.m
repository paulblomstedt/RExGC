function selectGenes(nTop)
% Selection of genes
% Author(s): Paul Blomstedt, Ritabrata Dutta

tabFilesRaw=dir('Data/TAB/Raw/*.tab');                      % files of gene expression measurement data
tabFilesAnalysed=dir('Data/TAB/Analysed/*.tab');            % files of experiment analytics (t-statistics and p-values)
load('Data/MAT/G44_EFO_combined.mat', 'expNames');          % list of accession numbers for experiments selected for the analysis


% Extract accession numbers from files of measurement data
expRaw = cell(length(tabFilesRaw),1);
for k = 1:length(tabFilesRaw)
    a = strsplit(tabFilesRaw(k).name,'-');
    expRaw(k) = {[char(a(2)) '-' char(a(3)) '-' char(a(4))]};
end
[~, indImportRaw] = intersect(expRaw,expNames);

% Create list of genes shared by all experiments
disp(['Extracting list of genes from experiment no. ' num2str(1)])
myData = importdata(['Data/TAB/Raw/' tabFilesRaw(indImportRaw(1)).name]);
geneListAll = myData.textdata(2:end,1);                                                     % list of genes shared by all experiments
for k = 2:length(indImportRaw)
    disp(['Extracting list of genes from experiment no. ' num2str(k)])
    myData = importdata(['Data/TAB/Raw/' tabFilesRaw(indImportRaw(k)).name]);
    geneListAll = intersect(myData.textdata(2:end,1),geneListAll,'sorted');                 % update genelist by intersecting through all experiments
end


% Extract accession numbers from files of experiment analytics
expAnalysed = cell(length(tabFilesAnalysed),1);
for k = 1:length(tabFilesAnalysed)
    a = strsplit(tabFilesAnalysed(k).name,'-');
    expAnalysed(k) = {[char(a(2)) '-' char(a(3)) '-' char(a(4))]};
end
[~, indImportAnalysed] = intersect(expAnalysed,expNames);

% For each experiment, choose the nTop first genes from the default ordering given by
% Expression Atlas, then take the union of these genes over all experiments.
nFiles = length(indImportAnalysed);
geneList = [];
for k = 1:nFiles
    disp(['Choosing the top ' num2str(nTop) ' genes from experiment no. ' num2str(k)])
    tStat = importdata(['Data/TAB/Analysed/' tabFilesAnalysed(indImportAnalysed(k)).name]);
    a = tStat.textdata(3:end,3);
    [~, I] = unique(a);
    u = a(sort(I));
    geneList = union(u(1:nTop),geneList);
end


% Create final list of genes and final data set
geneList = intersect(geneListAll,geneList);
dataProcessed = cell(1,length(expNames));
fracNaN = zeros(1,length(expNames));
for k=1:length(expNames)
    disp(['Creating processed dataset for experiment no. ' num2str(k)])
    myData = importdata(['Data/TAB/Raw/' tabFilesRaw(indImportRaw(k)).name]);
    [~,indSelect] = intersect(myData.textdata(2:end,1),geneList,'sorted');                      % find row indices for genes in the list
    dataProcessed{k} = myData.data(indSelect,:);
    fracNaN(k) = sum(sum(isnan(dataProcessed{k})))/numel(dataProcessed{k});                     % calculate fraction of missing values
    if fracNaN(k)>0
        [row,col] = find(isnan(dataProcessed{k}));
        rowmeans = nanmean(dataProcessed{k},2);
        for i = 1:length(row)
            for j = 1:length(col)
                dataProcessed{k}(row(i),col(j)) = rowmeans(row(i))...                           % impute if missing values exist
                    + randn * nanstd(dataProcessed{k}(row(i),:))*0.01;                          % perturb to make values unique
            end
        end
    end
end


save(['Data/MAT/44_top_' num2str(nTop) '.mat'],'dataProcessed','geneList','geneListAll','fracNaN');
fprintf(['\nSaved results in file Data/MAT/44_top_' num2str(nTop) '.mat\n\n'])

