function DEpValCorrelation
% Similarity matrix for DE approach
% Author(s): Paul Blomstedt, Ritabrata Dutta

% Read in experiment names 
load('Data/MAT/G44_EFO_combined.mat', 'expNames');                          % accession numbers of experiments used in the analysis
tabFilesAnalysed=dir('Data/TAB/Analysed/*.tab');                            % files of experiment analytics (t-statistics and p-values)
load('Data/MAT/44_top_5.mat','geneList')                  % list of genes to be used in the analysis

% Extract accession numbers from files of experiment analytics
expAnalysed = cell(length(tabFilesAnalysed),1);
for k = 1:length(tabFilesAnalysed)
    a = strsplit(tabFilesAnalysed(k).name,'-');
    expAnalysed(k) = {[char(a(2)) '-' char(a(3)) '-' char(a(4))]};
end
[~, indImportAnalysed] = intersect(expAnalysed,expNames);

% Extract analytics data and form differential expression profiles.
nExp = length(indImportAnalysed);
pVals = cell(2,nExp);
for k = 1:nExp
    disp(['Processing analytics data for experiment no. ' num2str(k)])
    indFiles = find(strcmp(expNames(k),expAnalysed)); 
    gene = [];
    pval = []; 
    for i = 1:length(indFiles)
        myData = importdata(['Data/TAB/Analysed/' tabFilesAnalysed(indFiles(i)).name]);
        a = myData.textdata(3:end,3);
        [~, I] = unique(a);
        geneTmp = a(sort(I));
        pvalTmp = myData.data(sort(I),2);
        gene = [gene;geneTmp];
        pval = [pval;pvalTmp];
    end
    [~,I] = unique(gene);
    gene = gene(sort(I));
    pval = pval(sort(I));
    pVals{1,k} = gene;
    pVals{2,k} = pval;
end

% Calculate correlations between profiles and form similarity matrix.
fprintf('\nCalculating correlations and forming the similarity matrix.\n\n')

% As an option to a pre-loaded list, use all available genes:
% geneList = intersect(pVals{1,1},pVals{1,2});
% for k = 1:nExp
%     geneList = intersect(geneList,pVals{1,k});
% end

pvalMatrix = zeros(length(geneList),nExp);
for k = 1:nExp
    [~,indPval] = intersect(pVals{1,k},geneList);
    pvalMatrix(:,k) = pVals{2,k}(indPval);
end

% D_spear = 1-squareform(pdist(pvalmatrix','spearman')); % optionally, use Spearman's correlation
D_corr = 1-squareform(pdist(pvalMatrix','correlation'));

save('Data/MAT/DE.mat','D_corr');
fprintf('\nSaved results in file Data/MAT/DE.mat\n\n')

