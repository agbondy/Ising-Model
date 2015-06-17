function matchclim(figs,varargin)
    % puts all heatmaps in the input list of figure handles into a common colormap,
    % bounded by the most extreme values across the entire set of heatmaps
    j=1;
    method='';
    clim=[];
    while j<=length(varargin)
        if strncmpi(varargin{j},'usedata',4)
            method='tight';
        elseif strncmpi(varargin{j},'clim',4)
            j=j+1;
            clim=varargin{j};
        end
        j=j+1;
    end  
    if isempty(clim)
        switch method
            case 'tight'
                cl=[];
                for f=1:length(figs)
                   imhs=imhandles(figs(f));
                   imdata = cat(1,imhs.CData);
                   cl=[cl ; imdata(:)];
                end               
            otherwise
                cl=[];
                for f=1:length(figs)
                    kids=get(figs(f),'children');
                    for k=1:length(kids)
                        if isprop(kids(k),'clim')
                            cl = [cl get(kids(k),'clim')];                
                        end
                    end
                end
        end
        clim=minmax(cl);
    end
    for f=1:length(figs)
        figure(figs(f));
        kids=get(gcf,'children');
        for k=1:length(kids)
            if isprop(kids(k),'clim')
                set(kids(k),'clim',clim);
            end
        end
    end
end
