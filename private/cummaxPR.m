function [mTprCorr,mPrecCorr]  = cummaxPR(mTpr,mPrec,mBaseline) 
% cummax for precision-recall curve
% Author: Paul Blosmtedt

if mTpr(1)~= 0
    mTprCorr = [0;mTpr];
    mPrecCorr = [0;mPrec];
end

mPrecCorr = mPrecCorr(end:-1:1);
if mPrecCorr(1)< mBaseline
    mPrecCorr(1) = mBaseline;
end
for i=2:length(mPrecCorr)
    if mPrecCorr(i)<mPrecCorr(i-1)
        mPrecCorr(i) = mPrecCorr(i-1);
    end
end
mPrecCorr = mPrecCorr(end:-1:1);

