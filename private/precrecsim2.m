function [mPrec, mTpr, mFpr, mThresh, mBaseline] = precrecsim2(score, sim)
% Wrapper function for prec_rec.
% Author: Paul Blomstedt, Ritabrata Dutta

nQueries = size(score,2);
prec = zeros(size(score));
tpr = ones(size(score));
fpr = ones(size(score));
thresh = zeros(size(score));
baseline = zeros(1,nQueries);

for i=1:nQueries
   [prec_tmp, tpr_tmp, fpr_tmp, thresh_tmp] = prec_rec(score(:,i), sim(:,i),'numThresh',nQueries);
   prec(1:length(prec_tmp),i) = prec_tmp;
   tpr(1:length(prec_tmp),i) = tpr_tmp;
   fpr(1:length(prec_tmp),i) = fpr_tmp;
   thresh(1:length(prec_tmp),i) = thresh_tmp;
   baseline(i) = sum(sim(:,i))/nQueries;
end

mPrec = mean(prec,2);
mTpr = mean(tpr,2);
mFpr = mean(fpr,2);
mThresh = mean(thresh,2);
mBaseline = mean(baseline);
