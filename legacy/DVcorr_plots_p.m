classdef DVcorr_plots_p < handle
methods(Hidden=true)
%% MAGR
%% ELIPSE
    function obj=plot_ellipse_p(obj)
        % PLOT ELLIPSE FOR BIVARIATE DECISION VARIABLE
        hFit=plotEllipse([obj.mu1 obj.mu2],[1 obj.rho; obj.rho 1],obj.CI,[],2,obj.colorEllipse); hold on
        axis(4.*[-1 1 -1 1]); axis square;
        plot([obj.cr1 obj.cr1],ylim,'k');
        plot(xlim,[obj.cr2 obj.cr2],'k');
        set(gca,'ytick',get(gca,'xtick'));

        % PLOT DECISION VARIALBE DATA (IF AVAILABLE)
        if min(obj.DATA(:)) ~= 0 || max(obj.DATA(:)) ~= 1
            plot(obj.DATA(:,1),obj.DATA(:,2),'ko','markerface','w','markersize',8);
        end
        Axis.format();
    end
    function obj=label_ellipse_p(obj)
        x='Decision Variable: Pass 1';
        y='Decision Variable: Pass 2';
        titl=['\rho=' num2str(obj.rho,'%.2f') ', \mu_1=' num2str(obj.mu1,'%.2f') ', \mu_2=' num2str(obj.mu2,'%.2f') ', c_1=' num2str(obj.cr2,'%.2f') ', c_2=' num2str(obj.cr2,'%.2f')];
        Axis.format(x,y,titl);
    end
    function obj=get_colors(obj)
        C=obj.nCmp;
        if ischar(obj.colormap)
            cmap=eval([obj.colormap ';']);
        else
            cmap=obj.colormap;
        end
        ind=round(linspace(1,size(cmap,1),C));
        obj.colors=cmap(ind,:);
    end
    function obj=get_color(obj,c)
        obj.colorEllipse=obj.colors(c,:);
    end
    function obj=label_cmp(obj)
        cmp=round(unique(obj.cmp),2);
        cmp=sprintf('%.2f',cmp);
        title(['Cmp ' cmp ]);
    end
%% RESPONSE
    function obj=plot_response_count_p(obj)
        % PLOT RESPONSE FREQUENCY
        S=obj.rspAgr;
        bar([S.Pnn S.Pnp S.Ppn S.Ppp],1,'w','linewidth',obj.LineWidth);
        xlim([0.25 4.75]);
        axis square;
        set(gca,'xtick',[1 2 3 4]);
        set(gca,'XTickLabel',{'--','-+','+-','++',});
        Axis.format();
    end
    function obj=label_response_count_p(obj)
        x='Joint Response Type';
        y='% Responses';
        titl='';
        Axis.format(x,y,titl);
    end
%% RHO PLOT
    function plot_rho_bin_p(obj,RHO,nBins)
        if ~exist('RHO','var')
            RHO=obj.RHO;
        end
        if ~exist('nBins','var') || isempty(nBins)
            X=Hist.bin_widths_FD(RHO);
        else
            X=nBins;
        end
        [counts,X]=hist(RHO,X);
        bar(X,counts,1,'w','LineWidth',obj.LineWidth);
    end
    function format_rho_bin_p(obj)
        Axis.format('Between pass correlation','Count');
    end

    function plot_rho_scatter_p(obj,RHO)
        if ~exist('RHO','var')
            RHO=obj.RHO;
        end
        x=repmat(unique(obj.stdX),size(RHO,1),1);
        plot(x,RHO,'ko','LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
    end
    function format_rho_scatter_p(obj)
        ylim([.01 3]);
        set(gca,'yscale','log');
        ytickformat('%.2f');
        yticks([.3 3]);
    end
    function xlabl=xlabel_rho_scatter_p(obj)
        if ~isempty(obj.Xname)
            Xname=obj.Xname;
        else
            Xname='X';
        end


        if ~isempty(obj.Xunits)
            xlabl=[' (' obj.Xunits ')'];
        else
            xlabl=[Xname];
        end
    end

    function ylabl=ylabel_rho_scatter_p(obj)
        ylabl='Between pass correlation';
    end

%% RATIO
    function plot_ratio_bin_p(obj,r,nBins)
        if ~exist('r','var')
            r=(1-obj.RHO)./(obj.RHO);
        end
        if ~exist('nBins','var')
            nBins=[];
        end
        %X=Hist.bin_widths_FD(log(r),1);
        %X=unique([X -X]);
        %[counts,ctrs]=hist(log(r),X);
        loghist(r,nBins,'LineWidth',obj.LineWidth,'EdgeColor',[0 0 0]);
        %set(gca,'Xtick',-1:1);
        %set(gca,'Xticklabel',10.^get(gca,'Xtick'));
        %xlim([-1,1]);
    end
    function format_ratio_bin_p (obj)
        Axis.format('\sigma^{2}_{Int}/\sigma^{2}_{Ext}','Count');
    end

    function plot_rho_scatter2_p(obj,sym,stagger)
        if ~exist('sym','var') || isempty(sym)
            sym='o';
        end
        if ~exist('r','var')
            r=obj.RHO;
        end
        if ~exist('stagger','var') || isempty(stagger)
            stagger=0;
        end
        x=repmat(unique(obj.stdX),size(r,1),1)+stagger;
        plot(x,r,['k' sym],'LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
    end

    function plot_ratio_scatter_p(obj,sym,stagger)
        if ~exist('sym','var') || isempty(sym)
            sym='o';
        end
        if ~exist('r','var')
            r=sqrt((1-obj.RHO)./(obj.RHO));
        end
        if ~exist('stagger','var') || isempty(stagger)
            stagger=0;
        end
        x=repmat(unique(obj.stdX),size(r,1),1)+stagger;
        plot(x,r,['k' sym],'LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
    end
    function format_ratio_scatter_p(obj)
        if ~isempty(obj.Xunits)
            xtitl=[obj.Xname ' (' obj.Xunits ')'];
        else
            xtitl=[obj.Xname];
        end
        Axis.format(xtitl,'\sigma^{2}_{Int}/\sigma^{2}_{Ext}');
        set(gca,'yscale','log');

        vals=[.1 .3  1  3 10];
        yticks(vals);
        yticklabels(vals);

        x=unique(obj.stdX);
        xticks(x(1));
        xticklabels(x);
    end
%% ABS
    function plot_abs_scatter_thresh_p(obj,sym)
        if ~exist('sym','var') || isempty(sym)
            sym='o';
        end
        I=sqrt(obj.TvarI);
        E=sqrt(obj.TvarE);
        T=sqrt(I^2+E^2);
        t=obj.T;
        x=unique(obj.stdX);

        plot(x,I,['b' sym],'LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        plot(x,E,['r' sym],'LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        %plot(x,T,['m' sym],'LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        %plot(x,t,'yd','LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        %plot(x,I./E,'ko','LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        %set(gca,'yscale','log')
        %legend('\sigma_I','\sigma_E','\sigma_T','\sigma_I/\sigma_E','location','northwest')
        legend('\sigma_i','\sigma_e','location','northwest');
    end
    function plot_abs_scatter_p(obj,color,bSplit)
        if ~exist('color','var') || isempty(color)
            color='k';
        end
        if ~exist('bSplit','var') || isempty(bSplit)
            bSplit=1;
        end
        i=obj.varI;
        e=obj.varE;
        if bSplit
            x=unique(obj.cmpX);
        else
            x=repmat(unique(obj.stdX),size(i,1),1);
        end
        plot(x,i,'ko','LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        plot(x,e,'rs','LineWidth',obj.LineWidth,'markerFaceColor','w','MarkerSize',10);
        if bSplit
            for j=1:length(x)
                plot([x(j),x(j)],[i(j),e(j)],'Color',color);
            end
        end
    end
    function format_abs_scatter_p(obj)
        if ~isempty(obj.Xunits)
            xtitl=[obj.Xname ' (' obj.Xunits ')'];
        else
            xtitl=[obj.Xname];
        end
        Axis.format(xtitl,'\hat{\sigma}^{2}_{Int} \hat{\sigma}^{2}_{Ext}');

        % XXX
        %set(gca,'yscale','log');
        %vals=[.1 .3  1  3 10];
        %yticks(vals);
        %yticklabels(vals);

        x=unique(obj.stdX);
        xticks(x(1));
        xticklabels(x);
    end

end
end
