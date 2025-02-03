classdef Magr_plot < handle
methods
    function plotMagr(obj,opts)
        if nargin > 2
            opts=[];
        end
        obj.plot_magr(opts);
        obj.format_magr(opts);
    end
    function plot_magr(obj,opts,fitType)
        if nargin < 2
            opts=[];
            bPlotBinom=true;
        elseif isstruct(opts) && isfield(opts,'bPlotBinom')
            bPlotBinom=opts.bPlotBinom;
        else
            bPlotBinom=true;
        end
        if nargin < 3
            if isstruct(opts) && isfield(opts,'fitType')
                fitType=opts.fitType;
            else
                fitType='func';
            end
        end
        % NULL MODEL
        if bPlotBinom
            obj.plot_magr_binom_parab(opts); hold on
        end

        switch fitType
        case 'func'
            obj.plot_magr_fit_func(opts); hold on
        case 'parab'
            obj.plot_magr_fit_parab(opts); hold on
        case 'sim'
            obj.plot_magr_sim_boot(opts); hold on
        end

        % DATA
        obj.plot_magr_data(opts);
        obj.format_magr(opts);
    end
    function obj=plot_magr_binom_parab(obj,opts)
        if isempty(obj.PBino)
            obj.get_agree_binom();
        end
        P=obj.PBino;

        opts=Plot.parse(opts,'FillColor',[0 0 0], ...
                             'FaceAlpha',.3, ...
                             'EdgeColor','none', ...
                             'LineStyle','--', ...
                             'LineWidth',1.5, ...
                             'LineColor',[0 0 0], ...
                             'CI',68);
        if opts.CI==68
            U=P.PA68U;
            L=P.PA68L;
        elseif opts.CI==95
            U=P.PA95U;
            L=P.PA95L;
        end

        [xb,yb,Xb,Yb]=obj.get_agree_interp(P.PC,P.PA,U,L);

        % plot binom CI
        fill(Xb,Yb, opts.FillColor, ...
                    'FaceAlpha',opts.FaceAlpha,  ...
                    'EdgeColor',opts.EdgeColor);
        hold on;
        plot(xb,yb,'Color',opts.LineColor', ...
                   'LineStyle','--', ...
                   'LineWidth',opts.LineWidth);
        hold off;
    end
    function obj=plot_magr_fit_func(obj,opts)
        if isempty(obj.PFit)
            obj.get_agree_fit();
        end
        if isempty(obj.PSim)
            obj.get_agree_sim();
        end

        FillColor=[];
        LineColor=[];

        opts=Plot.parse(opts,'FillColor',[0 0 0], ...
                             'FaceAlpha',.3, ...
                             'EdgeColor','none', ...
                             'LineStyle','-', ...
                             'LineWidth',1.5, ...
                             'LineColor',[0 0 0], ...
                             'CI',68);

        P=obj.PSim;
        if opts.CI==68
            U=P.PA68U;
            L=P.PA68L;
        elseif opts.CI==95
            U=P.PA95U;
            L=P.PA95L;
        end

        % FILL
        [xb,yb,Xb,Yb]=obj.get_agree_interp(P.PC,P.PA,L,U);
        fill(Xb,Yb,opts.FillColor, ...
                   'FaceAlpha',opts.FaceAlpha, ...
                   'EdgeColor',opts.EdgeColor);
        hold on;
        % LINE
        [PC,idx]=sort(obj.PFit.PC);
        PA=obj.PFit.PA(idx);
        plot(P.PC,P.PA, ...
                 'LineStyle',opts.LineStyle, ...
                 'Color',opts.LineColor, ...
                 'LineWidth',opts.LineWidth);
        hold off;

    end
    function obj=plot_magr_sim_boot(obj,opts)
        %if isempty(obj.pA)
        %end

        P=PSim;
        opts=Plot.parse(opts,'FillColor',[0 0 0], ...
                             'FaceAlpha',.3, ...
                             'EdgeColor','none', ...
                             'LineStyle','-', ...
                             'LineWidth',1.5, ...
                             'LineColor',[0 0 0], ...
                             'CI',68 ...
                             );
        if opts.CI==68
            U=P.PA68U;
            L=P.PA68L;
        elseif opts.CI==95
            U=P.PA95U;
            L=P.PA95L;
        end

        [xb,yb,Xb,Yb]=obj.get_agree_interp(P.PC,P.PA,L,U);
        fill(Xb,Yb,opts.FillColor, ...
                   'FaceAlpha',opts.FaceAlpha, ...
                   'EdgeColor',opts.EdgeColor);
        hold on;
        [PC,idx]=sort(obj.PSim.PC);
        PA=obj.PSim.PA(idx);
        plot(100*PC,100*PA, ...
                 'LineStyle',opts.LineStyle, ...
                 'Color',opts.LineColor, ...
                 'LineWidth',opts.LineWidth);
        hold off;
    end
    function obj=plot_magr_fit_parab(obj,opts)
        if isempty(obj.PCfit)
            obj.get_fit_parab();
        end


        k=1000;
        X=linspace(-.5,.5,k);
        x=X+.5;
        % XXX
        y=polyval(obj.PCfit,X);
        yu=polyval(obj.PCfitU,X);
        yl=polyval(obj.PCfitL,X);

        X=[x fliplr(x)];
        Y=[yl fliplr(yu)];

        opts=Plot.parse(opts,'FillColor',[0 0 0], ...
                             'FaceAlpha',.3, ...
                             'EdgeColor','none', ...
                             'LineStyle','-', ...
                             'LineWidth',1.5, ...
                             'LineColor',[0 0 0]);

        patch(X,Y, ...
             opts.FillColor, ...
            'FaceAlpha',opts.FaceAlpha, ...
            'EdgeColor',opts.EdgeColor); hold on;
        plot(x,y, ...
                 'LineStyle',opts.LineStyle, ...
                 'Color',opts.LineColor, ...
                 'LineWidth',opts.LineWidth);
        hold off;
    end
    function obj=plot_magr_data(obj,opts)
        if isempty(obj.PEmp)
            obj.get_agree_emp();
        end
        opts=Plot.parse(opts, ...
                             'Marker','o',...
                             'MarkerFaceColor',[1 1 1], ...
                             'MarkerEdgeColor',[0 0 0], ...
                             'MarkerSize',15, ...
                             'LineStyle','none', ...
                             'LineWidth',1.5, ...
                             'LineColor',[0 0 0], ...
                             'FontSize',15);

        plot(obj.PEmp.PC,obj.PEmp.PA, ...
                         'LineStyle', 'none', ...
                         'LineWidth', opts.LineWidth, ...
                         'Marker',    opts.Marker, ...
                         'MarkerSize',opts.MarkerSize, ...
                         'MarkerEdgeColor',opts.MarkerEdgeColor, ...
                         'MarkerFaceColor',opts.MarkerFaceColor);

        if numel(unique(obj.FIT.RHO))==1
            txt=['\rho =' num2str(obj.FIT.RHO(1),'%0.2f')];
            text(.3, .95,txt,'FontSize',opts.FontSize);
        end
        hold off;
    end
    function obj=format_magr(obj,opts)
        if nargin < 2 || isempty(opts)
            opts=struct();
        end
        if ~isfield(opts,'bTitle')
            opts.bTitle=true;
        end
        xlbl='% Comparison Chosen';
        ylbl='% Agreement';
        if opts.bTitle
            titl=['stdX = ' num2str(round(obj.stdXunq(1),2,'significant'))];
        else
            titl=[];
        end
        Axis.format(xlbl,ylbl,titl);
        axis square;
        ylim([0.4 1]);
        set(gca,'XTick',(0:0.5:1));
    end
end
end
