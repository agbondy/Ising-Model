function makeCorrelationImage(varargin)
    warning('off','MATLAB:hg:willberemoved');
    warning('off','MATLAB:warn_r14_structure_assignment');
    p=inputParser;
    p.KeepUnmatched=true;
    p.addOptional('ax',gca,@(x)validateattributes(x,{'matlab.graphics.axis.Axes'},{'scalar'}));
    p.addParamValue('labels',{},@(x)validateattributes(x,{'cell','char'},{'nonempty'}));
    p.addParamValue('title',{},@(x)validateattributes(x,{'cell','char'},{'nonempty'}));
    p.parse(varargin{:});
    xl=get(p.Results.ax,'xlim');
    yl=get(p.Results.ax,'ylim');
    if isprop(gca,'xticklabelrotation')
        set(p.Results.ax,'xtick',ceil(xl(1)):floor(xl(2)),'ytick',ceil(yl(1)):floor(yl(2)),...
            'xticklabel',p.Results.labels,'yticklabel',p.Results.labels,'xaxislocation','top',...
            'xticklabelrotation',35,'yticklabelrotation',0);
    else
        set(p.Results.ax,'xtick',ceil(xl(1)):floor(xl(2)),'ytick',ceil(yl(1)):floor(yl(2)),...
            'xticklabel',p.Results.labels,'yticklabel',p.Results.labels,'xaxislocation','top');        
    end
    h=colorbar;
    if isfield(get(h),'Label')
        h.Label.String='correlation';
    end
    axis square;
    title(p.Results.title);    
end
