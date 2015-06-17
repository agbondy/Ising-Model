function makeCorrelationImage(varargin)
    warning('off','MATLAB:hg:willberemoved');
    p=inputParser;
    p.KeepUnmatched=true;
    p.addOptional('ax',gca,@(x)validateattributes(x,{'matlab.graphics.axis.Axes'},{'scalar'}));
    p.addParameter('labels',{},@(x)validateattributes(x,{'cell','char'},{'nonempty'}));
    p.addParameter('title',{},@(x)validateattributes(x,{'cell','char'},{'nonempty'}));
    p.parse(varargin{:});
    xl=get(p.Results.ax,'xlim');
    yl=get(p.Results.ax,'ylim');
    set(p.Results.ax,'xtick',ceil(xl(1)):floor(xl(2)),'ytick',ceil(yl(1)):floor(yl(2)),...
        'xticklabels',p.Results.labels,'yticklabels',p.Results.labels,'xaxislocation','top',...
        'xticklabelrotation',35,'yticklabelrotation',0);
    h=colorbar;
    h.Label.String='correlation';
    axis square;
    title(p.Results.title);    
end
