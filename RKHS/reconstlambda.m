function [lamvec, shatvec] = reconstlambda(k,xhat,ktvec,tt)
%RECONSTLAMBDA Summary of this function goes here
%   Detailed explanation goes here
    assert(sum(size(xhat)-size(ktvec)).^2==0);
    shatvec = k(ktvec(:),tt(:).').'*xhat;
    lamvec = exp(shatvec);
end