function rexVisualization(EFOtype,algorithm,nTop)
% Visualization of experiment retrieval using gene clustering.
% Author(s): Paul Blomstedt, Ritabrata Dutta
%
% Inputs:
% EFOtype: 'cell_type', 'disease' or 'organism_part'
% algorithm: cell describing the choice of retrieval method
% nTop: identifier for different datasets based on the choice of genes

%%%%%


% Load the saved results and plot Prec-Rec curve
load(['Data/MAT/G44_' EFOtype '.mat'],'G','indKeep');
G(logical(eye(size(G)))) = 0; % set diagonal to 0

figure
[~,~,~,~,mBaseline] = precrecsim2(G, G);                            % baseline for Prec-Rec curve
plot([0 1], [mBaseline mBaseline],'k-.')

hold all
for i = 1:length(algorithm)
    if any(strcmp(algorithm{i},{'PPM','PPMkm','k-means'}))
        load(['Data/MAT/' algorithm{i} '_NID_top_' num2str(nTop) '.mat'],'D_nid');
        dM = 1-D_nid;
        dM = dM(indKeep,indKeep);
        dM(logical(eye(size(dM)))) = 0;
        [mPrec, mTpr] = precrecsim2(dM,G);
        [mTpr,mPrec]  = cummaxPR(mTpr,mPrec,mBaseline);
        plot(mTpr,mPrec)
    elseif strcmp(algorithm{i},'LogML')
        load(['Data/MAT/' algorithm{i} '.mat'],'D_lml');
        dM = D_lml;
        dM = dM(indKeep,indKeep);
        dM(logical(eye(size(dM)))) = 0;
        [mPrec, mTpr] = precrecsim2(dM,G);
        [mTpr,mPrec]  = cummaxPR(mTpr,mPrec,mBaseline);
        plot(mTpr,mPrec)
    else
        load(['Data/MAT/' algorithm{i} '.mat'],'D_corr');
        dM = D_corr;
        dM = dM(indKeep,indKeep);
        dM(logical(eye(size(dM)))) = 0;
        [mPrec, mTpr] = precrecsim2(dM,G);
        [mTpr,mPrec]  = cummaxPR(mTpr,mPrec,mBaseline);
        plot(mTpr,mPrec)
    end
end
hold off

lim = get(gca,'ylim');
ylim([max(0,lim(1)-0.02) min(1,lim(2)+0.02)])

xlabel('Recall','FontSize',16)
ylabel('Precision','FontSize',16)
set(gca,'FontSize',16)






