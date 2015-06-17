function x = shuffle(x)
    % like randperm but output will not return ANY indices that are the
    % same as the input, i.e. removes all simultaneity in time series
    if ~ismatrix(x) && length(x)>1
        mssg(0,'Only works on non-singleton vectors.');
        return;
    end
    success=0;
    len=length(x);
    inds=1:len;
    while ~success
        newinds=randperm(len);
        if sum(newinds==inds)==0
            success=1;
        end
    end
    x=x(newinds);
end