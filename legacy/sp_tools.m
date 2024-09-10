classdef sp_tools < handle
properties
    EXPS
    sp
    varargin

    subjSel
    passesSel
    stdSel

    nPssAtOnce
    nSubj
    nPss
    nStd
    nCmp

    nSubjPlot
    nPssPlot
    nStdPlot

    I
    subj
    nInd
    bCmp
    bComb
    bSame
    MODE

    Opts
    spOpts
    mainSpOpts
    ITER
    iter
    continueflag
end
methods
    function obj=sp_tools(PAR,bCmp,bSame,bComb,Opts,varargin)

        obj.nPssAtOnce=PAR.nPssAtOnce;
        obj.nPss=PAR.nPss;
        obj.nSubj=PAR.nSubj;
        obj.nStd=PAR.nStd;
        obj.nCmp=PAR.EXPS.nCmp;

        obj.EXPS=PAR.EXPS;

        obj.bCmp=bCmp;
        obj.bComb=bComb;
        obj.bSame=bSame;
        obj.Opts=Opts;

        obj.get_IND(PAR);

        obj.parse_varargin(varargin);
        obj.get_sp_opts();
        obj.parse_opts();
        obj.get_sp();

    end
    function obj=parse_varargin(obj,vin)
        if isempty(vin)
            return
        end

        if length(vin) > 0 && ~isempty(vin{1})
            obj.subjSel=vin{1};
        end
        if length(vin) > 1 && ~isempty(vin{2})
            obj.passesSel=vin{2};
            if iscell(obj.passesSel) && ~isempty(obj.passesSel)
                obj.passesSel=cellfun(@transpose_fun,obj.passesSel,'UniformOutput',false);
                obj.passesSel=vertcat(obj.passesSel{:});
            elseif isnumeric(obj.passesSel)
                if size(obj.passesSel,1)==obj.nPssAtOnce && size(obj.passesSel,2)~=obj.nPssAtOnce
                    obj.passesSel=transpose(obj.passesSel);
                elseif size(obj.passesSel,2)~=obj.nPssAtOnce
                    error('incorrect dimensions for nPssAtOnce & passes');
                end
            end
        else
            obj.passesSel=[1,2];
        end
        if length(vin) > 2 && ~isempty(vin{3})
            obj.stdSel=vin{3};
        end

        function out=transpose_fun(in)
            if size(in,1)>1 & size(in,2)==1
                out=transpose(in);
            else
                out=in;
            end
        end
    end
    function obj=get_sp_opts(obj)
        Opts=struct();
        Opts=obj.get_sp_size(Opts);
        Opts=obj.get_sp_ticks(Opts);
        Opts=obj.get_sp_titles(Opts);
        obj.spOpts=Opts;
    end

    function Opts=get_sp_titles(obj,Opts)
        if obj.bCmp && nflds(obj.EXPS) > 0
            Opts.ctitl=obj.EXPS.get_rc_title('cmp');
            Opts.rtitl=obj.EXPS.get_rc_title('std');
            %Opts.sub=obj.EXPS.get_rc_title('subjs', obj.subjSel);
        elseif obj.bCmp
            Opts.rtitl=strcat('Std ',   strobj.split(Num.toStr([1:obj.nStd]),','));
            Opts.sub=num2str(obj.subj);
        elseif nflds(obj.EXPS) > 0
            Opts.xtitl=obj.EXPS.get_xtitl();
            Opts.rtitl=obj.EXPS.get_rc_title('imgDim','pass');
            if obj.bSame
                Opts.ctitl=[];
            else
                Opts.ctitl=obj.EXPS.get_rc_title('subj');
            end
            Opts.rtitl=join(nchoosek(Opts.rtitl,obj.nPssAtOnce));
        else
            Opts.xtitl=obj.cell{1}.xlabel_rho_scatter_p();
            Opts.rtitl=strcat('Passes ', strobj.split(Num.toStr(obj.ITER.PASSES),';'));
            if bSame
                Opts.citl=[];
            else
                Opts.ctitl=strcat('Subj ',   strobj.split(Num.toStr([1:obj.nSubj]),','));
            end
        end
    end
    function Opts=get_sp_size(obj,Opts)
        if obj.bSame
            nsubj=1;
        elseif ~isempty(obj.subjSel) && isnumeric(obj.subjSel)
            nsubj=length(obj.subjSel);
        else
            nsubj=numel(unique(obj.ITER.IND(:,2)));
        end
        obj.nSubjPlot=nsubj;

        if ~isempty(obj.passesSel) && iscell(obj.passesSel)
            npss=length(obj.subjSel)
        elseif ~isempty(obj.passesSel) && isnumeric(obj.passesSel)
            npss=1;
        else
            npss=obj.nPss/obj.nPssAtOnce;
        end
        obj.nPssPlot=npss;

        if ~isempty(obj.stdSel) && isnumeric(obj.stdSel)
            nstd=length(obj.stdSel);
        else
            nstd=obj.nStd;
        end
        obj.nStdPlot=nstd;

        if obj.bCmp
            Opts.sz=[nsubj nstd];
        else
            Opts.sz=[npss nsubj];
        end
    end
    function Opts=get_sp_ticks(obj,Opts)
        Opts.xticks=unique(obj.EXPS.stdXunqAll(:));
        Opts.xtickAngle=40;
    end
    function obj=parse_opts(obj)
        fldsP=properties('subPlots');
        fldsM={'sz','xtitl','ytitl','titl','rtitl','ctitl'};

        % OVERWITE SPOPTS WITH EXISTING OPTS
        fldsS=fieldnames(obj.spOpts);
        flds=fieldnames(obj.Opts);
        for i = 1:length(flds)
            fld=flds{i};
            if ismember(fld,fldsS) && ~isempty(obj.Opts.(fld))
                obj.spOpts.(fld)=obj.Opts.(fld);
            end
        end

        % POPULATE ALL
        ALL=obj.Opts;
        obj.Opts=struct();

        flds=fieldnames(obj.spOpts);
        for i = 1:length(flds)
            fld=flds{i};
            ALL.(fld)=obj.spOpts.(fld);
        end
        obj.spOpts=struct();
        obj.mainSpOpts=struct();

        % DISTRIBUTE ALL
        flds=fieldnames(ALL);
        for i = 1:length(flds)
            fld=flds{i};
            if ismember(fld,fldsM)
                obj.mainSpOpts.(fld)=ALL.(fld);
            elseif ismember(fld,fldsP)
                obj.spOpts.(fld)=ALL.(fld);
            else
                obj.Opts.(fld)=ALL.(fld);
            end
        end

        % POPULATE MISSING MAIN
        fldsS=fieldnames(obj.spOpts);
        fldsm=fieldnames(obj.mainSpOpts);
        for i=1:length(fldsM)
            fld=fldsM{i};
            if ~ismember(fld,fldsm)
                obj.mainSpOpts.(fld)=[];
            end
        end
        if isfield(obj.Opts,'syms')
            obj.ITER.syms=obj.Opts.syms;
            obj.Opts=rmfield(obj.Opts,syms);
        end
        if isfield(obj.Opts,'bStagger') && obj.Opts.bStagger
            obj.ITER.staggers=obj.get_staggers(obj.ITER.X);
            obj.Opts=rmfield(obj.Opts,'bStagger');
        end
        if isfield(obj.Opts,'colors') && ~isempty(obj.Opts.colors)
            obj.ITER.colors=obj.Opts.colors;
            obj.Opts=rmfield(obj.Opts,'colors');
        end
        if ~isfield(obj.Opts,'bFlipRC')
            obj.Opts.bFlipRC=0;
        end

    end
    function obj=get_IND(obj,PAR)
        obj.ITER=struct();

        obj.ITER.PASSES=PAR.get_passes();
        if obj.bCmp==2
            obj.ITER.IND=PAR.get_ind_cmp();
        else
            obj.ITER.IND=PAR.get_ind_ind();
        end
        obj.nInd=size(obj.ITER.IND,1);

        obj.ITER.seen=[];
        if obj.bCmp
            % XXX
        else
            obj.ITER.X=unique(obj.EXPS.stdXunqAll(:));
        end
        if obj.bSame & obj.bComb
            obj.MODE=3;
        elseif obj.bComb
            obj.MODE=2;
        elseif obj.bSame
            obj.MODE=1;
        else
            obj.MODE=0;
        end
        obj.get_I();
    end
    function obj=get_I(obj)
        % std  = 1
        % subj = 2
        % cmp  = 3
        % ind  = end
        ONES=true(size(obj.ITER.IND,1),1);
        if isempty(obj.subjSel)
            SUBJ=ONES;
        else
            SUBJ=ismember(obj.ITER.IND(:,2),obj.subjSel);
        end

        if isempty(obj.passesSel)
            PASS=ONES;
        else
            % XXX
            PASS
        end

        if isempty(obj.stdSel)
            STD=ONES;
        else
            STD=ismember(obj.ITER.IND(:,1),obj.stdSel);
        end
        obj.I=SUBJ & PASS & STD;

    end
    function obj=get_iter(obj,i,R)
        obj.iter=struct();

        obj.get_inds(i);
        obj.get_continueflag();
        if obj.continueflag==1; return; end
        obj.get_seen();
        obj.sp_select(i);

        % SYMS
        if isfield(obj.ITER,'syms') && ~isempty(obj.ITER.syms);
            obj.iter.sym=obj.ITER.syms(obj.iter.subj);
        end

        % STAGGERS
        if isfield(obj.ITER,'staggers') && ~isempty(obj.ITER.staggers);
            obj.iter.stagger=obj.ITER.staggers(obj.iter.std, obj.iter.subj);
        end
        if isfield(obj.ITER,'colors') && ~isempty(obj.ITER.colors);
            obj.iter.color=obj.ITER.colors(obj.iter.subj);
        end

        % R
        if ~exist('R','var') || isempty(R)
            if ~isfield(obj.iter,'X') && isfield(obj.ITER,'X')
                obj.iter.X=obj.ITER.X;
            end
            return
        end
        bAll=obj.MODE==3;
        obj.get_R_iter_ind(i,R,bAll);

    end
    function obj=sp_select(obj,i)
        % bCMP
        %   nsubj nstd
        % ELSE
        %   npss nsubj
        % bSame -> column

        % HANDLE SAME
        if obj.bCmp
            r=obj.get_subj_plot_ind();
        else
            r=obj.get_pass_plot_ind();
        end

        if obj.bSame
            c=1;
        elseif obj.bCmp
            c=obj.get_std_plot_ind();
        else
            c=obj.get_subj_plot_ind();
        end


        % SELECT
        if obj.Opts.bFlipRC
            obj.sp.select(c,r);
        else
            obj.sp.select(r,c);
        end
    end
    function out=get_subj_plot_ind(obj)
        if ~isempty(obj.subjSel) && isnumeric(obj.subjSel)
            out=find(obj.subjSel==obj.iter.subj);
        else
            out=obj.iter.subj;
        end

    end
    function out=get_pass_plot_ind(obj)
        if ~isempty(obj.passesSel) && isnumeric(obj.passesSel)
            out=find(ismember(obj.passesSel, obj.iter.passes, 'rows'));
        else
            out=obj.iter.passes;
        end
    end
    function out=get_std_plot_ind(obj)
        % STD SEL
        if ~isempty(obj.stdSel) && isnumeric(obj.stdSel)
           out=find(obj.stdSel==obj.stdSel);
       else
           out=obj.iter.std;
        end
    end
    function obj=get_inds(obj,i)
        % std  = 1
        % subj = 2
        % cmp  = 3
        % ind  = end
        obj.iter=struct;
        obj.iter.IND=obj.ITER.IND(i,:);
        obj.iter.std=obj.ITER.IND(i,1);
        obj.iter.subj=obj.ITER.IND(i,2);
        obj.iter.ind=obj.ITER.IND(i,end);
        if obj.bCmp==2
            obj.iter.cmp=obj.ITER.IND(i,4);
        end
        obj.iter.passes=obj.ITER.PASSES(obj.ITER.IND(i,3),:);
    end
    function obj=get_continueflag(obj)
        % CORRECT SUBJ AND PASS
        obj.continueflag=0;

        % SUBJ SEL
        bA=~isempty(obj.subjSel) && isnumeric(obj.subjSel) && ~isequal(obj.iter.subj, obj.subjSel);
        if bA
            obj.continueflag=1;
        end

        bA=~isempty(obj.passesSel) && isnumeric(obj.passesSel) && ~ismember(obj.iter.passes, obj.passesSel,'rows');
        if bA
            obj.continueflag=1;
        end

        % STD SEL
        bA=~isempty(obj.stdSel) && isnumeric(obj.stdSel) && ~isequal(obj.iter.std, obj.stdSel);
        if bA
            obj.continueflag=1;
        end

        % SEEN
        obj.iter.bSeen=obj.bComb && ismember(obj.iter.subj,obj.ITER.seen);
        if obj.iter.bSeen
            obj.continueflag=1;
        end

    end
    function obj=get_seen(obj);

        % UPDATE
        obj.ITER.seen=[obj.ITER.seen obj.iter.subj];
    end
    function obj=get_R_iter_ind(obj,i,R,bAll)

        flds=fieldnames(R);
        for k = 1:length(flds)
            fld=flds{k};
            if bAll && length(R.(fld))==obj.EXPS.nStd && numel(R.(fld))==length(R.(fld));
                obj.iter.(fld)=R.(fld);
            elseif obj.MODE==2
                ind=obj.ITER.IND(:,2)==obj.iter.subj;
                obj.iter.(fld)=R.(fld)(ind);
            elseif size(R.(fld),3) > 1
                obj.iter.(fld)=R.(fld)(i,:,:);
            else
                obj.iter.(fld)=R.(fld)(i);
            end
        end
        if ~isfield(obj.iter,'X')
            obj.iter.X=obj.ITER.X;
        end
    end

    function obj=get_sp(obj)
        if obj.Opts.bFlipRC
            r=obj.mainSpOpts.rtitl;
            c=obj.mainSpOpts.ctitl;
            obj.mainSpOpts.rtitl=c;
            obj.mainSpOpts.ctitl=r;
            obj.mainSpOpts.sz=fliplr(obj.mainSpOpts.sz);
        end
        sz=obj.mainSpOpts.sz;
        xtitl=obj.mainSpOpts.xtitl;
        ytitl=obj.mainSpOpts.ytitl;
        titl=obj.mainSpOpts.titl;
        rtitl=obj.mainSpOpts.rtitl;
        ctitl=obj.mainSpOpts.ctitl;
        opts=obj.spOpts;

        obj.sp=SubPlots(sz,xtitl,ytitl,titl,rtitl,ctitl,opts);
    end
    function plot_interv(obj)
         Plot.interv(obj.iter.X, obj.iter.U, obj.iter.L); hold on;
    end
    function plot(obj)
        plot(obj.iter.X, obj.iter.m,'ko', 'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
    end

    function obj=finalize(obj)
        obj.sp.finalize();
    end
    function [m,U,L]=average_stat(obj,val,bAll,bRmNeg)
        if ~exist('bRmNeg','var') || isempty(bRmNeg)
            bRmNeg=0;
        end
        % std subj pssind
        % ind cmp
        if bRmNeg
            val(val<0)=nan;
            val=log10(val);
        end
        if size(val,1) > 1 & size(val,2) > 1
            M=nanmean(val,2);
        else
            M=val;
        end
        if bAll
            for i = 1:obj.EXPS.nStd
                ind=obj.ITER.IND(:,1)==i;
                m(i)=mean(M(ind));
                stdev(i)=nanstd(M(ind));
            end
        else
            m=M;
            stdev=nanstd(val,[],2);
        end

        U=m+stdev;
        L=m-stdev;
        if bRmNeg
            m=10.^(m);
            U=10.^U;
            L=10.^L;
        end
    end
    function staggers=get_staggers(obj,X)
        d=diff(X)/2;
        staggers=cell2mat(arrayfun(@(x) linspace(0,x,obj.nSubj)-x/2,d,UniformOutput,false));
        staggers=[staggers(1,:); staggers];
    end
end
end
