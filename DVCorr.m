classdef DVCorr < handle & DVcorr_plots & DVcorr_plots_p & Magr
properties
    DVFITTER

    nParabSim
    nTrlParabSim

    % plot params
    colorEllipse
    colormap
    CIcolor
    LineWidth
    Xname
    Xunits
    colors

    % fitt dv corr
    FIT
    BOOT

    % FROM CURVES
    % COMBINED PASSES
    varI
    varE
    varT

    bProg=true

    RCMPCHS
    TBL
    CMPX
    CMPXI
    CMPXO
    STDX
    stdXunq
    cmpXunq
    nCmp
    stddev1
    stddev2
end
methods
    function obj= DVCorr(stdX2,cmpX2,RcmpChs2,Opts,stddev1,stddev2)
    %function obj= DVCorr(stdX2,cmpX2,RcmpChs2,Opts)

        % XXX after dataTable
        bTest=nargin < 3;
        if nargin >=1 && ~isempty(stdX2)
            obj.STDX=stdX2;
        end
        if nargin >=2 && ~isempty(cmpX2)
            CMPX=cmpX2;
        end
        if nargin >= 3 && ~isempty(RcmpChs2);
            obj.RCMPCHS=RcmpChs2;
        end
        if nargin < 4 || isempty(Opts)
            Opts=struct();
        end
        if size(RcmpChs2,2) == 1
            error('Only one pass detected');
        end


        obj.DVFITTER=DVFitter(obj);
        obj.parse(Opts);

        if bTest
            obj.RCMPCHS = mvnrnd([0 0],[1 .5; .5 1],1000);
        end

        % RM_NANS
        %try
            ii=any(isnan(obj.RCMPCHS),2) | any(isnan(CMPX),2) | any(isnan(obj.STDX),2);
            if any(ii)
                obj.RCMPCHS(ii,:)=[];
                obj.STDX(ii,:)=[];
                CMPX(ii,:)=[];
            end
        %end
        obj.CMPXO=CMPX;
        [~,~,obj.CMPXI]=unique([CMPX obj.STDX],'rows');
        obj.CMPX=obj.CMPXO(:,1);

        obj.stdXunq=unique(obj.STDX,'rows');
        obj.cmpXunq=unique(obj.CMPX);
        obj.nCmp=numel(obj.cmpXunq);
    end
    function Opts=parse(obj,Opts)

        P=DVCorr.getP();
        [S,Opts]=Args.simpleLoose([],P,Opts);
        Args.applyIf(obj,S);

        Opts=Opts';
        Opts=Opts(:);
        Opts=obj.DVFITTER.parse(Opts);


    end
%- INIT
    function getStdDevs(obj)
        copts=struct('nBoot',1,'nBest',obj.DVFITTER.nBest);
        curve1=psyCurve(obj.STDX(:), obj.CMPXO(:), obj.RCMPCHS(:),copts);

        obj.stddev1=sqrt(curve1.tFit.^2 * curve1.DPcrt);
        obj.stddev2=obj.stddev1;
    end
    function getStdDevsBoth(obj)
        copts=struct('nBoot',1,'nBest',obj.DVFITTER.nBest);
        curve1=psyCurve(obj.STDX(:,1), obj.CMPXO(:,1), obj.RCMPCHS(:,1),copts);
        curve2=psyCurve(obj.STDX(:,2), obj.CMPXO(:,2), obj.RCMPCHS(:,2),copts);

        obj.stddev1=sqrt(curve1.tFit.^2 * curve1.DPcrt);
        obj.stddev2=sqrt(curve2.tFit.^2 * curve2.DPcrt);
    end
    function run(obj,bProg)
        if strcmp(obj.DVFITTER.modelType,'R') && isempty(obj.stddev1)
            obj.getStdDevs();
        elseif strcmp(obj.DVFITTER.modelType,'R*') && isempty(obj.stddev1)
            obj.getStdDevsBoth();
        end

        if nargin < 2; bProg=obj.bProg; end
        assignin('base','DVCORR',obj);
        obj.DVFITTER.run(bProg);
        obj.get_magr();
    end
    function out=get.FIT(obj)
        out=obj.DVFITTER.Fit;
    end
    function out=get.BOOT(obj)
        out=obj.DVFITTER.Boot;
    end
% XXX SELECT
    function obj=get_noise_sources(obj,PsyCurve)
        %C=obj.nCmp;
        %obj.varE=zeros(C,1);
        %obj.varI=zeros(C,1);
        %obj.varT=zeros(C,1);


        %for c = 1:C
        %    obj.select(c);
        %    dp=PsyCurve.DPfit(c);
        %    if obj.cmp==obj.std
        %        continue
        %    end
        %    [obj.varE(c),obj.varI(c),obj.varT(c)]=DVcorr.getNoiseSourcesInd(obj.cmp,obj.std,obj.rho,dp);
        %end


        dp=PsyCurve.DPcrt;
        T=PsyCurve.tFit;
        rho=obj.rho;
        [obj.varE,obj.varI,obj.varT]=DVCorr.getNoiseSources(T,rho,dp);
    end
end
methods(Static=true)
    function P=getP()
        P={ ...
            'colorEllipse', [0 0 0],'';
            'colormap',     'cmapGray2','';
            'CIcolor',      0.9.*[1 1 1],'';
            'LineWidth',    2,'';
            'Xname',        [],'';
            'Xunits',       [],'';

            'nParabSim',    1000,'';
            'nTrlParabSim', [],'';
        };
    end
    function [TvarE,TvarI,TvarT]=getNoiseSources(T,rho,dp)
        varT=T^2;
        varE=rho*TvarT;
        varI=TvarT-TvarE;


        %Rideal=mu1Fix
        %mu/sigma
        [~,~,~,~,DPfit,DPdta]=psyfitgengaussAll(stdX,cmpX,Rideal,stdXunq,[],[1],1.36,0);
        clear mu1Fix;
        %for s = 1:length(unique(abs(stdX)))
        disp('psyFitDecisionVariableCorrAll: WARNING! mu1Fix is klooged!!! BE CAREFUL!!!');
        mu1Fix(:,s) = DPfit(:,s).*sqrt(.43)./sqrt(2);
        obj.stddev2=sqrt(curve2.tFit.^2 * curve2.DPcrt);

    end
    function out=genData(stdX,cmpX,Mu,sigma,rho,cr,nTrlPerCmp)
        nPass=size(sigma,1);
        nCmp=size(cmpX,2);

        % RESIZE
        % Cmpx
        cmpXAll=repelem(cmpX(1,:)',nTrlPerCmp,1);
        if numel(stdX)==1
            stdXAll=repmat(stdX,size(cmpXAll));
        else
            TODO
        end

        if isempty(Mu)
            Mu=cmpX;
        end
        if size(Mu,1)==1
            Mu=repmat(Mu,nPass,1);
        end
        if size(cr,2)==1
            cr=repmat(cr,1,nCmp);
        end

        Rho=rho2Rho(rho);
        Sigma=rho2cov(rho,sigma);
        crit=(Mu-cr)';


        %  DAT
        DVDat=[];
        RDat=[];
        for i = 1:nCmp
            % DV DATA
            dvdat=mvnrnd(Mu(:,i),Sigma,nTrlPerCmp);
            DVDat=[DVDat; dvdat];

            % R DATA
            rdat =bsxfun(@gt,dvdat,crit(i,:));
            RDat=[RDat; rdat];
        end
        out=struct();
        out.RCmpChs=RDat;
        out.DV=DVDat;
        out.cmpX=cmpXAll;
        out.stdX=stdXAll;
        out.Sigma=Sigma;
        out.Rho=rho;

    end
end
end
