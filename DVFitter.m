classdef DVFitter < handle
% DVFit
%     RHO(ndim,ncmp)
% std dev
% combinecmp
properties
    DVFit
    PsyFit

%% DATA
    RCMPCHS
    STDX
    CMPX
    CMPXI

%% Options

    %: MODEL SPEC
    modelType
    extType
    bObs
    % AND/OR specify *Inds
    rhoFitInd % fit Inds
    muFitInd
    crFitInd
    stdFitInd
    %
    rhoFix
    muFix
    crFix
    stdFix
    %
    bCombineCmp % combine
    bCtrMu
    nSplit

    cls
    bNoFit

    %: Repeat OPTS
    nBest
    nBoot
    bootCIPrcnt

    bProg

%% CALC
    nCmp

    %: fmincon
    fminOpts
    mvnopts
    mvnconsts
    F00

    %: DIMENSIONS
    nPass
    nAry
    bBinary
    nQuad
    nRhoAll
    nAgree
    %
    patterns
    agreeInds

    %: Param Constants
    nRhoFit
    nMuFit
    nCrFit
    nStdFit
    %
    bRhoFit
    bMuFit
    bCrFit
    bStdFit
    %
    bRhoFitAny
    bMuFitAny
    bCrFitAny
    bStdFitAny
    bMuCrFitAny
    %
    rhoFitIndParam % fmincon input to expanded params
    muFitIndParam
    crFitIndParam
    stdFitIndParam
    %
    LB
    UB

    %: reepeat structs
    Fit

    bootSeed

end
properties(Hidden)
    bBoot
    bBest
    Parent
    optsFlds
    opts
    inopts


    % current fit
    boot
    best

    RCMPCHS_
    CMPX_
    STDX_
    CMPXI_
end
methods(Static)
    function obj=new(stdX,cmpX,RCmpChs,Opts)
        obj=DVFitter();
        obj.STDX=stdX;
        obj.CMPX=cmpX;
        obj.RCMPCHS=RCmpChs;
        [R,~,obj.CMPXI]=unique([cmpX stdX],'rows');

        obj.parse(Opts);
    end
    function P=getP()
        P={ ...
            'modelType',    'RMM','';
            ...
            'rhoFitInd',       [],'';
            'muFitInd',        [],'';
            'crFitInd',        [],'';
            ...
            'rhoFix',       [],'';
            'muFix',        [],'';
            'crFix',        [],'';
            'stdFix',        [],'';
            ...
            'bCombineCmp',  true,'';
            'nSplit',   0,'';
            'bCtrMu',  true,'';
            ...
            'bNoFit',       false,'';
            'cls',          'double','';
            'bObs',         false,'';
            ...
            'nBoot',        1000,'';
            'nBest',        100,'';
            ...
            % fminopts
            'minFuncType',  'fmincon','';
            'bParallel',    0,'';
            'bootSeed',     1,'';
            'bootCIPrcnt',  68.27,'';
            'bProg',        true,'';
        };
    end
    function out=getDefaults()
        P=DVFitter.getP();
        out=P(:,1:2)';
        out=struct(out{:});
    end
end
methods
%- GET/SET
    function out=get.bBoot(obj)
        out=obj.nBoot > 1;
    end
    function out=get.bBest(obj)
        out=obj.nBest > 0;
    end
    %
    function set.RCMPCHS(obj,val)
        if ~isempty(obj.Parent)
            obj.Parent.RCMPCHS=val;
        else
            obj.RCMPCHS_=val;
        end
    end
    function out=get.RCMPCHS(obj)
        if ~isempty(obj.Parent)
            out=obj.Parent.RCMPCHS;
        else
            out=obj.RCMPCHS_;
        end
    end
    %
    function set.CMPX(obj,val)
        if ~isempty(obj.Parent)
            obj.Parent.CMPX=val;
        else
            obj.CMPX_=val;
        end
    end
    function out=get.CMPX(obj)
        if ~isempty(obj.Parent)
            out=obj.Parent.CMPX;
        else
            out=obj.CMPX_;
        end
    end
    %
    function set.CMPXI(obj,val)
        if ~isempty(obj.Parent)
            obj.Parent.CMPXI=val;
        else
            obj.CMPXI_=val;
        end
    end
    function out=get.CMPXI(obj)
        if ~isempty(obj.Parent)
            out=obj.Parent.CMPXI;
        else
            out=obj.CMPXI_;
        end
    end
    %
    function set.STDX(obj,val)
        if ~isempty(obj.Parent)
            obj.Parent.STDX=val;
        else
            obj.STDX_=val;
        end
    end
    function out=get.STDX(obj)
        if ~isempty(obj.Parent)
            out=obj.Parent.STDX;
        else
            out=obj.STDX_;
        end
    end
    %
%- CONSTRUCT
    function obj=DVFitter(parent)
        if nargin < 1
            return
        end
        obj.Parent=parent;
    end
%- PARSE
    function Opts=parse(obj,Opts)
        P=DVFitter.getP();
        if iscell(Opts)
            [S,Opts]=Args.simpleLoose([],P,Opts{:});
        else
            [S,Opts]=Args.simpleLoose([],P,Opts);
        end
        [~,opts]=Args.applyIf(obj,S);

        obj.inopts=opts;

        % fminopts

    end
    function get_fminOpts(obj,minfunctype)
        if nargin < 2
            minfuncType='fmincon';
        end
        alg='sqp';
        bSparse=false;
        %Display='final';
        Display='off';
        FDType='central';
        bUseParallel=false;
        maxIter=500;
        Args={...
            'Algorithm',alg;
            'UseParallel',bUseParallel;
            'Display',Display;
            'MaxIterations',maxIter;
            'StepTolerance',1e-12;
            'MaxFunctionEvaluations',1000;
            'OptimalityTolerance',1e-10;
            'FiniteDifferenceType',FDType;
            'ConstraintTolerance',1e-10;

        };
        switch alg
        case 'active-set'
            Args2={...
                'FuncitonTolerance',10^-10;
            };
        case 'interior-point'
            Args2={...
                'SubproblemAlgorithm','factorization';
                'HonorBounds',true;
            };
        case 'sqp'
            Args2={...
                'ObjectiveLimit',-1e20;
                %'ScaleProblem',  false;
                'ScaleProblem',  'none';
            };
        otherwise
            Args2={};
        end

        Args=[Args; Args2]';
        F=optimoptions(minfunctype,Args{:});

        obj.fminOpts=F;
    end
%- RUN
    function cont(obj,bProg)
        if nargin < 2
            bProg=obj.bProg;
        end

        % MAIN FIT
        if isempty(obj.best)
            obj.run_best(bProg);
        end
        if obj.nBoot > 1
            obj.run_boot(bProg,true);
        end
    end
    function run(obj,bProg,bContinue)
        if nargin < 2
            bProg=obj.bProg;
        end
        if nargin < 3 || isempty(bContinue)
            bContinue=false;
        end
        if ~bContinue
            obj.init();
        end

        % best
        obj.run_best(bProg);

        % boot
        if obj.nBoot > 1
            obj.run_boot(bProg,false);
        end

    end
%- INIT
    function init(obj)
        obj.optsFlds=obj.getFitOptsFlds();

        % dims
        DVFitter.get_dims(obj);

        % params
        obj.init_params();

        % mvnopts
        obj.mvnopts=rMvnCdf.getOpts(obj.nPass);
        if obj.nPass > 3
            obj.mvnopts.bCheckBounds=false;
        end

        % mvnconstants
        obj.mvnconsts=rMvnCdf.getConstants(obj.nPass,obj.cls,obj.mvnopts);

        % F0
        obj.F00=DVFitter.getF0(obj);

        obj.get_fminOpts(obj.inopts.minFuncType);
        obj.getOpts();

        % PACK
        obj.init_pack;
    end
    function init_params(obj)

        if obj.bCombineCmp
            tol=1e-2;
        else
            tol=0.05;
        end

        if isempty(obj.modelType)
            mdl='';
        else
            mdl=upper(strrep(obj.modelType,'*',''));
        end
        if contains(mdl,{'E','F'})
            obj.extType=mdl;
            ext=mdl;
            switch mdl
            case {'E','F','E0','F0'}
                mdl='RR';
            case 'EL'
            end
        else
            obj.extType='';
            ext=mdl;
        end

        if obj.bCtrMu
            muDef=-1;
        else
            muDef=0;
        end
        % fld, n, bound, default
        flds={...
             'rho', obj.nRhoAll, [0 1],0;
             'mu' , obj.nPass,4,muDef;
             'cr' , obj.nPass,3,0;
             'std', obj.nPass,nan,1;
        };

        LB=cell(3,1);
        UB=cell(4,1);
        count=0;
        for i = 1:size(flds,1)
            % flds
            fld=flds{i,1};
            f=fld(1);
            F=upper(f);
            Fld=[F fld(2:end)];

            n=flds{i,2};
            b=flds{i,3};
            defFix=flds{i,4};

            fInd=[fld 'FitInd'];
            bFFit=['b' Fld 'Fit'];
            bFFitAny=['b' Fld 'FitAny'];
            nFFit=['n' Fld 'Fit'];
            fIndParam=[fld 'FitIndParam'];

            fFix=[fld 'Fix'];

            % fit indeces
            if ~isempty(obj.(fInd))
                nIn=numel(obj.(fInd));
                if nIn~=n
                    error(['Number of elements in field ' fld 'Ind must be ' num2str(n) '.']);
                elseif Vec.isRow(obj.(fInd))
                    obj.(fInd)=obj.(fInd)';
                end
            elseif contains(mdl,[F F])
                obj.(fInd)=(1:n)';
            elseif contains(mdl,F)
                obj.(fInd)=ones(n,1);
            else
                obj.(fInd)=zeros(n,1);
            end

            % bFFit
            obj.(bFFit)=obj.(fInd) > 0;
            obj.(bFFitAny)=any(obj.(bFFit));

            % fIndParam
            %obj.(nFFit)=sum(obj.(bFFit));
            obj.(nFFit)=sum(unique(obj.(fInd))>0);
            if (obj.(bFFitAny))
                obj.(fIndParam)=count+(1:obj.(nFFit));
                count=count+obj.(nFFit);
            end

            % fix values
            fixval=ones(size(obj.(fInd)))*defFix;
            fixval(obj.(fInd)>0)=nan;
            fixset=obj.(fFix);
            if any(~isnan(fixset)) && any(obj.(fInd)==0)
                fixval(obj.(fInd)==0)=fixset;
            end
            obj.(fFix)=fixval;

            % BOUNDS
            nB=obj.(nFFit);
            if nB==0
                continue
            end
            if numel(b)>1
                LB{i}=ones(nB, 1)* b(1) + tol;
                UB{i}=ones(nB, 1)* b(2) - tol;
            else
                LB{i}=ones(nB, 1)*-b + tol;
                UB{i}=ones(nB, 1)* b - tol;
            end
        end
        obj.LB=vertcat(LB{:})';
        obj.UB=vertcat(UB{:})';

        obj.bMuCrFitAny=obj.bMuFitAny || obj.bCrFitAny || obj.bStdFitAny;
    end
%- PACK
    function init_pack(obj,flds,pflds)
        % rho               = (nRhoAll,ncmp)
        % mu,cr,std         = (ndim,ncmp)
        % NegLL,            = (1,ncmp)
        % NegLLAll,         = (1,1)
        % %
        % Ntrl,PC           = (1,ncmp)
        % NtrlAll           = (1,1)
        % PA                = (nAgree,nCmp)
        % P                 = (nQuad,nCmp)

        [flds,pflds]=obj.get_pack_flds();
        F=struct();
        if obj.bBoot
            b=struct();
            nCmp=numel(unique(obj.CMPX));
        end

        for i = 1:length(flds)
            fld=flds{i};
            F.(fld)=struct();
            if strcmp(fld,'P')
                continue
                for i = 1:length(pflds)
                    F.P.(pfld)=struct('best',[]);
                    %F.P.(pfld)=struct('best',[],'ci',[],'std',[],'mean',[],);
                end
            else
                F.(fld)=struct('best',[],'ci',[],'std',[],'mean',[]);

                if obj.bBoot

                    switch fld
                    case 'NegLLAll'
                        n=1;
                        m=1;
                    case 'NegLL'
                        n=1;
                        m=nCmp;
                    case 'RHO'
                        n=obj.nRhoAll;
                        m=nCmp;
                    case 'STD'
                        n=obj.nPass;
                        m=1;
                    otherwise
                        n=obj.nPass;
                        m=nCmp;
                    end
                    b.(fld)=zeros(n,m,obj.nBoot);
                end
            end
        end
        obj.Fit=F;
        if obj.bBoot
            obj.boot=b;
        end

    end
    function best_pack(obj)
        [flds,pflds]=obj.get_pack_flds();
        F=obj.Fit;
        for i =1:length(flds)
            fld=flds{i};
            if strcmp(fld,'P')
                for j = 1:length(pflds)
                    pfld=pflds{j};

                    F.P.(pfld).best=obj.DVFit.P.(pfld);
                end
            else
                F.(fld).best=obj.DVFit.(fld);
            end
        end
        obj.Fit=F;
    end
    function boot_pack(obj,b)
        [flds,pflds]=obj.get_pack_flds();
        for i =1:length(flds)
            fld=flds{i};
            if strcmp(fld,'P')
                continue
                %for j = 1:length(pflds)
                %    pfld=pflds{i};
                %    obj.boot.P.(pfld)(:,:,b)=obj.DVFit.P.(pfld);
                %end
            else
                obj.boot.(fld)(:,:,b)=obj.DVFit.(fld);
            end
        end
    end
    function boot_pack_complete(obj)
        [flds,pflds]=obj.get_pack_flds();

        F=obj.Fit;
        d=3; % boot dim
        for i = 1:length(flds)
            fld=flds{i};

            if strcmp(fld,'P')
                for j = 1:length(pflds)
                    pfld=(pflds{j});

                    %val=obj.boot.P.(pfld)( ~isnan(any(obj.boot.P.(pflds{j}),d)),:);
                    %F.P.(pfld).ci=confInt(val,obj.bootCIPrcnt,d);
                    %F.P.(pfld}).std=std(val,[],d);
                    %F.P.(pfld).mean=mean(val,d);
                end
            else

                val=obj.boot.(fld);
                F.(fld).ci  =confInt(val, obj.bootCIPrcnt,d);
                F.(fld).std =std(val,[],d);
                F.(fld).mean=mean(val,d);
            end
        end
        obj.Fit=F;
    end
    function [flds,pflds]=get_pack_flds(obj)
        flds={'RHO', 'MU', 'CR','STD','NegLL','NegLLAll','P'};
        pflds={'nTrl','nTrlAll','PC','PA','P','N'};
    end
    function getOpts(obj)
        flds=obj.optsFlds;
        obj.opts=struct();
        for i = 1:length(flds)
            obj.opts.(flds{i})=obj.(flds{i});
        end
    end
%- FIT ROUTINES
    function run_test(obj)
        obj.fminOpts.Display='final';
        obj.DVFit=DVFit(obj.opts);
        obj.DVFit.applyData(obj.RCMPCHS, obj.CMPX, obj.CMPXI, obj.STDX);
        obj.safe_run();
    end
    function run_boot(obj,bProg,bContinue)
        if nargin < 2 || isempty(bProg)
            bProg=true;
        end
        if nargin < 3 || isempty(bContinue)
            bContinue=any(obj.boot.NegLLAll(:) > 0);
        end
        bContinue=bContinue && ~isempty(obj.boot);

        % progress
        if ~bContinue
            strt=1;
        else
            strt=find(obj.boot.NegLLAll(1,1,:)==0,1,'first');
        end
        if bProg
            p=Pr(obj.nBoot-strt+1,1,'Bootstrapping');
        end

        % OPTS
        obj.getOpts();

        % RND
        rng(obj.bootSeed,'twister');
        sds=randi(2^32-1,obj.nBoot,1);
        N=size(obj.RCMPCHS,1);

        obj.DVFit=DVFit(obj.opts);
        for i=strt:obj.nBoot
            if bProg
                p.u();
            end

            rng(sds(i),'twister');
            ii = randi(N,N,1);

            obj.DVFit.applyData(obj.RCMPCHS(ii,:),obj.CMPX(ii), obj.CMPXI(ii,:), obj.STDX(ii));
            obj.safe_run;

            obj.boot_pack(i);
        end
        if bProg
            p.c();
        end
        obj.boot_pack_complete();
    end
    function run_best(obj,bProg,bContinue)
        if nargin < 2 || isempty(bProg)
            bProg=true;
        end
        if nargin < 3 || isempty(bContinue)
            bContinue=~isempty(obj.best);
        end
        bProg=bProg && obj.nBest > 1;

        % PROGRESS
        if ~bContinue
            strt=1;
            obj.best.NegLL=zeros(obj.nBest,1);
            best=inf(1);
        else
            strt=find(obj.best.NegLL==0,1,'first');
            best=min(obj.best.NegLL(obj.best.NegLL>0));
        end
        if bProg && obj.nBest > 1
            p=Pr(obj.nBest-strt+1,1,'Finding Best');
        end

        % OPTS
        obj.getOpts();

        % INIT DVFIT
        obj.DVFit=DVFit(obj.opts);
        obj.DVFit.applyData(obj.RCMPCHS,obj.CMPX, obj.CMPXI, obj.STDX);

        for i = strt:obj.nBest
            if bProg && obj.nBest > 1
                p.u();
            end

            % FIT
            obj.safe_run();

            % SELECT IF BETTER
            ind=obj.DVFit.NegLLAll < best;
            obj.best.NegLL(i)=obj.DVFit.NegLLAll;
            if ~any(ind)
                continue
            end
            obj.best_pack();

        end
        if bProg && obj.nBest > 1
            p.c();
        end
    end
    function safe_run(obj)
        count=0;
        while true
            try
                obj.DVFit.run();
                break
            catch ME
                count=count+1;
                if count > 5
                    rethrow(ME);
                end
            end
        end
    end
%- HELPERS
end
methods(Static)
    function flds=getFitOptsFlds()
        fldsP=props(DVFitter);
        flds=props(DVFit);
        flds=flds(ismember(flds,fldsP));
    end
    function varargout=get_dims(obj)
        bObj=isa(obj,'DVFitter');
        if bObj
            R=obj.RCMPCHS;
        else
            R=obj;
            obj=struct();
        end

        nPass=size(R,2)/max([1,obj.nSplit]);
        nAry=numel(unique(R(:)));
        nQuad=nAry^nPass; %number of orthants
        if size(R,3)==1
            nRhoAll=ceil((nQuad-nPass)/2);
        else
            nRhoAll=ceil((nQuad-nPass)/2)+nPass;
        end
        bBinary=nAry==2;

        % All patterns and agreeInds
        [patterns,agreeInds]=Magr.getRespPats(R,nPass);
        nAgree=max(agreeInds);


        if nargout <= 1
            obj.nPass=nPass ;
            obj.nAry=nAry;
            obj.nQuad=nQuad;
            obj.nRhoAll=nRhoAll;
            obj.bBinary=bBinary;
            obj.patterns=patterns;
            obj.agreeInds=agreeInds;
            obj.nAgree=nAgree;
            varargout={obj};
        else
            TODO
        end
    end
    function F0=getF0(obj)
        z1=zeros(obj.nPass,1);
        z2=zeros(obj.nPass);
        [L0,U0]=getCRLU0(obj);
        L=L0;
        U=U0;
        L(~isinf(L))=0;
        U(~isinf(U))=0;
        e=eye(obj.nPass);
        t=ones(obj.nPass);
        bDiag=size(obj.RCMPCHS,3)~=1;
        if ~bDiag;
            indRhoU=~e & triu(t);
        else
            indRhoU=triu(t)==1;
        end

        %[cr1 cr2]
        %[cr1 Inf]
        %[Inf cr2]
        if obj.nPass==2
            CRS=[0; 1; 1];
        else
            CRS=[];
        end

        F0=struct('CRL',L0,...
                  'CRU',U0,...
                  'CRL0',L,...
                  'CRU0',U,...
                  ...
                  'CRB',U0,...
                  'CRS',CRS,...
                  ...
                  'Mu',z1,...
                  ...
                  'indRhoU', indRhoU,...
                  'Rho',z2,...
                  'RhoU',z2,...
                  'bRhoD',bDiag,...
                  'teye',e==1,...
                  'eye',e...
        );

        function [L,U]=getCRLU0(obj)
            [~,Y]=ndgrid(1:obj.nQuad,1:obj.nPass);
            pat=repmat({1:obj.nAry},1,obj.nPass);
            pat=Set.distribute(pat{:});
            L=pat;
            U=pat;
            iL=L==1;
            iU=U~=1;
            L(iL)=-inf;
            U(iU)=inf;
            for i = 1:obj.nPass
                L(~iL & Y==i)=i;
                U(~iU & Y==i)=i;
            end
            [~,ind]=sort(sum(~isinf(L),2));
            L=L(ind,:);
            U=U(ind,:);
        end
    end
end
end
