function hh = mysgtitle(varargin)
    if  verLessThan('matlab', '9.5')
        title = varargin{1};
        ax = setAxis(gcf,1);
        if not(isempty(ax.Title.String))
            if iscell(ax.Title.String)
                ax.Title.String = [{title} ax.Title.String];
            else
                ax.Title.String = [{title} {ax.Title.String}];
            end
        else
            ax.Title.String = {title};
        end
    else
        sgtitle(varargin(:));
    end
end
