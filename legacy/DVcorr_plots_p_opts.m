classdef DVcorr_plots_p_opts < handle
methods
    function Opts=get_ellipse_opts(obj,Opts,sub)
        Opts.xtitl='1st Pass Decision Variable';
        Opts.ytitl='2nd Pass Decision Variable';

        str=obj.get_sub_str(sub);
        Opts.titl=['Observer ' sub ' DV correlation'];
        Opts.ax='square';
        Opts=Struct.combinePref(Opts,Opts);
    end
    function Opts=get_counts_opts(obj,Opts,sub)
        nCmp=obj.get_nCmp();
        Opts.sz=[obj.nStd nCmp];
        Opts.ytitl='Proportion Chosen';
        Opts.xtitl='Joint Reobj.sponse pattern';

        % Opts.xtickLabels=pm(xtickInd); XXX
        % Opts.xticks=1:size(Opts.xtickLabels,1); XX
        Opts.marginDE=[6 6 11 6];

        str=obj.get_sub_str(sub);
        Opts.titl=['Observer ' sub ' Pass Agreement'];
        Optsylim=[.4 1]; % NOTE
        Opts.ax='square';

        Opts=Struct.combinePref(Opts,Opts);
    end
    function Opts=get_magr_opts(obj,Opts,sub)
        Opts.ytitl='Proprtion Cmp. Chosen';
        Opts.xtitl='Proportion Agreement';
        str=obj.get_sub_str(sub);
        Opts.titl=['Observer' str ' pass Agreement'];
        Opts.xticks=[0 .25 .5 .75 1];
        Opts.yticks=[.5 .75 1];

        Opts.ax='square';
        Opts=Struct.combinePref(Opts,Opts);
    end
    function Opts=get_bin_ratio_opts(obj,opts,sub)
        str=obj.get_sub_str(sub);
        Opts.xtitl='Importance Ratio';
        Opts.ytitl='Count';
        Opts.titl=['Observer ' sub ' Ratio Counts'];
        Opts.bLog=1;
        Opts.ax='square';


        Opts=Struct.combinePref(opts,Opts);
    end
    function  Opts=get_bin_rho_opts(obj,opts,sub)
        str=obj.get_sub_str(sub);
        Opts.xtitl='Between pass Correlation';
        Opts.ytitl='Count';
        Opts.titl=['Observer ' sub ' Correlation Counts'];
        Opts.bLog=0;
        Opts.ax='square';

        Opts=Struct.combinePref(opts,Opts);
    end
    function Opts=get_thresh_opts(obj,Opts)
        Opts.sz=[1,1];
        Opts.xtitl='Pass 1 Thresholds (arcmin)';
        Opts.ytitl='Pass 2 Thresholds (arcmin)';
        Opts.titl='Discrimination Thresholds';
        Opts.yscale='log'
        Opts.xscale='log'
        Opts.ax='square';
        Opts=Struct.combinePref(opts,Opts);
    end
    function Opts=get_ratio_scatter_opts(obj,opts)
        Opts.ytitl='Importance Ratio';
        Opts.titl='Internal noise to stimulus variability importance';
        Opts.yscale='log';
        Opts.ax='square';

        Opts=Struct.combinePref(opts,Opts);
    end
    function Opts=get_rho_scatter_opts(obj,opts)
        Opts.ytitl='Between pass correlation';
        Opts.titl='DV correlation';
        Opts.ax='square';
        Opts=Struct.combinePref(opts,Opts);
    end
    function Opts=scatter_abs_opts(obj,Opts)
        Opts.titl='Stimulus variability and internal noise estimates';

        Opts.ax='square';
        Opts=Struct.combinePref(Opts,Opts);

    end
    function str=get_sub_str(obj,sub)
        if iscell(sub) && numel(sub)==1
            str=[' ' sub{1}];
        elseif iscell(sub)
            str=[' ' join(subj,'-')];
        else
            str='';
        end
    end
end
end
