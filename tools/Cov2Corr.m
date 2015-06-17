function ExpCorr = Cov2Corr(ExpCovariance,ExpVar)
    if nargin<2
        ExpVar = diag(ExpCovariance);
    end
    ExpVar=ExpVar(:);
    ExpCorr = ExpCovariance ./ sqrt(ExpVar * ExpVar');
    ExpCorr(isnan(ExpCovariance))=NaN;
end