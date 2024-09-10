classdef DVFit < handle & DVFit_objFuns & DVExt
% NOTE RUN IN DVFit_fitRoutineS
properties
%% fitted paramas
        RHO
        MU
        MUSTD
        CR
        STD
        NegLL
        NegLLAll
        P

%% DATA FROM PARENT
        RCmpChs
        Rlbl
        stdX
        cmpX
        cmpXInd

        bCtrMu
        bExtConstr=0
        nSplit

%% RECALCULTED FROM DATA
        stdXunq
        cmpXunq
        cmpXIndunq
        nCmp


%% OPTS FROM PARENT
        bNoFit
        extType
        bObs

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
        bCombineCmp % NEGLL FROM SINGLE CMP OR TOGETHER

        %: fmincon
        fminOpts
        mvnopts
        mvnconsts
        F00
        In0

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
        bMuCrFitAny
        bStdFitAny
        %
        rhoFitIndParam % fmincon input to expanded params
        muFitIndParam
        crFitIndParam
        stdFitIndParam
        %
        LB
        UB

        %: fit iter
        raw
        nllFun
        %: constraints
        lb
        ub
        nlcon
        A
        b
        param0

        rho
        mu
        cr
        st
        mustd
        %
        p
        negLL
        negLLAll

        out0

end
properties(Hidden)
    %: sel
    C
    rcmpchs
    rlbl
    rcmps
    cmp
    cmpR
    cmpI
    ncmp
    cmpXunqI
    cmpXunqC

    muFixSet


end
methods
%- Init
    function obj=DVFit(Opts, varargin)
        if nargin < 1
            return
        end

        % OPTS

        flds=fieldnames(Opts);
        for i = 1:length(flds)
            obj.(flds{i})=Opts.(flds{i});
        end

        obj.mvnopts.TolFun=1e-12;
        % DATA
        if nargin < 1
            obj.applyData(varargin{:});
        end
    end
    function applyData(obj,RCmpChs,cmpX,cmpXInd,stdX)

        %% DATA
        % handle inputs
        obj.RCmpChs=RCmpChs;
        obj.cmpX=cmpX;
        obj.cmpXInd=cmpXInd;
        obj.stdX=stdX;

        obj.stdXunq=unique(obj.stdX);
        obj.cmpXunq=unique(obj.cmpX);
        obj.cmpXIndunq=unique(obj.cmpXInd);
        obj.nCmp=numel(obj.cmpXunq);

        % label according to patterns

        nCross=size(obj.RCmpChs,3);
        if obj.nSplit
            for j = 1:obj.nCross
            for i = 1:obj.nSplit
                inds=(1:obj.nPass)+(i-1)*obj.nPass;
                [~, obj.Rlbl{i}(:,j)]=ismember(obj.RCmpChs(:,inds,j),obj.patterns,'rows');
            end
            end
        else

            for j = 1:nCross
                [~, obj.Rlbl(:,j)]=ismember(obj.RCmpChs(:,:,j),obj.patterns,'rows');
            end
        end
        % XXX packageing


    end
%- Fit routines
    function run(obj)
        obj.init_package();
        if obj.bCombineCmp
            obj.fit_combine();
        else
            obj.fit_each();
        end
    end
    function fit_main(obj,bCombine)
        obj.get_bounds();
        obj.get_nlcon();

        % In0
        obj.center_mu();
        obj.get_In0(bCombine);

        % func
        obj.get_out0(); % OUT
        obj.get_objective();

        obj.get_param0(); % IN

        % Minimize
        if obj.bNoFit
            param=param0;
        else
            [obj.raw,~] = fmincon(obj.nllFun,obj.param0,[],[],[],[],obj.lb,obj.ub,obj.nlcon,obj.fminOpts);
        end
    end
    function fit_each(obj);
        C=obj.nCmp;
        for c = 1:C
            IND=all(obj.cmpXInd==obj.cmpXIndunq(c),2);
            obj.sel_cmp(IND);

            obj.fit_main(false);

            obj.package_each(c);
        end
    end
    function fit_combine(obj)
        obj.sel_cmp();
        obj.fit_main(true);
        obj.package_combine();
    end
%- Select
    function sel_cmp(obj,cmpInd)
        if nargin < 2 || isempty(cmpInd)
            cmpInd=1:size(obj.cmpX,1);
            bStable=true;
        else
            bStable=false;
            uarg={}; % XXX stable here too?
        end
        % SEL cmp
        obj.cmp=obj.cmpX(cmpInd,:);
        obj.cmpI=obj.cmpXInd(cmpInd,:);
        obj.cmpR=repmat(obj.cmp,1,obj.nPass);

        if bStable
            obj.cmpXunqI=unique(obj.cmpI,'stable');
            obj.cmpXunqC=unique(obj.cmp,'stable');
            [obj.cmpXunqI,ind]=sort(obj.cmpXunqI);
            obj.cmpXunqC=obj.cmpXunqC(ind);
        else
            obj.cmpXunqI=unique(obj.cmpI);
            obj.cmpXunqC=unique(obj.cmp);
        end
        obj.ncmp=numel(obj.cmpXunqC);

        % SEL R
        obj.rcmpchs=obj.RCmpChs(cmpInd,:,:);
        if obj.nSplit
            for i = 1:obj.nSplit
                obj.rlbl{i}=obj.Rlbl{i}(cmpInd,:);
            end
        else
            obj.rlbl=obj.Rlbl(cmpInd,:);
        end

        nCross=size(obj.RCmpChs,3);
        % get rcmps
        if obj.nSplit
            obj.rcmps=cell(obj.ncmp,obj.nSplit);
            for i=1:obj.nSplit
            for c = 1:obj.ncmp
                ind=all(obj.cmpI==obj.cmpXunqI(c),2);
                obj.rcmps{c,i}=obj.rlbl{i}(ind,:);
            end
            end
        else
            obj.rcmps=cell(obj.ncmp,1);
            for c = 1:obj.ncmp
                ind=all(obj.cmpI==obj.cmpXunqI(c),2);
                obj.rcmps{c}=obj.rlbl(ind,:);
            end
        end
    end
%- Pack
    function init_package(obj)
        if obj.bCombineCmp
            C=1;
        else
            C=obj.nCmp;
        end
        obj.NegLL    = zeros(1,C);
        obj.NegLLAll = zeros(1,1);

        obj.RHO      = zeros(obj.nRhoAll,C);
        obj.MU       = zeros(obj.nPass,C);
        obj.CR       = zeros(obj.nPass,C);
        obj.STD       = zeros(obj.nPass,C);
    end
    function package_each(obj,c);
        [obj.negLLAll,obj.negLL,obj.p]=obj.nllFun(obj.raw);

        if c == 1
            obj.P    = obj.p;
        else
            obj.P(c) = obj.p;
        end
        obj.NegLL(1,c)    = obj.negLL;
        obj.NegLLAll(1,1) = sum(obj.NegLL,2);
        obj.RHO(:,c)      = obj.rho;
        obj.MU(:,c)       = obj.mu;
        obj.CR(:,c)       = obj.cr;
        obj.STD(:,c)      = obj.cr;

    end
    function package_combine(obj)
        [obj.negLLAll,obj.negLL,obj.p]=obj.nllFun(obj.raw);
        if contains(obj.extType,'E')
            F=obj.get_params_ext(obj.raw);
            raw=F.Rho;
        else
            raw =obj.raw;
        end

        z=zeros(obj.nPass,obj.ncmp);
        obj.MU=z;
        obj.CR=z;
        obj.RHO=zeros(obj.nRhoAll,obj.ncmp);
        for c = 1:obj.ncmp
            obj.get_params(raw,c);
            obj.MU(:,c)=obj.mu;
            obj.CR(:,c)=obj.cr;
            if obj.nSplit
                % HERE
            else
                obj.RHO(:,c)=obj.rho;
            end
        end
        obj.STD      = obj.st;

        obj.rho      = obj.RHO;
        obj.mu       = obj.MU;
        obj.cr       = obj.CR;
        obj.st       = obj.STD;

        obj.P        = obj.p;
        obj.NegLL    = obj.negLL;
        obj.NegLLAll = obj.negLLAll;
    end
%- Params
    function get_params(obj,param,c)
        rho=obj.rhoFix;
        if obj.bRhoFitAny
            p=param(obj.rhoFitIndParam);
            %rho(obj.bRhoFit)=p(obj.rhoFitInd); NOTE CHANGED
            rho(obj.bRhoFit)=p(obj.rhoFitInd(obj.bRhoFit));
        end

        mu =obj.muFixSet; % NOTE
        if obj.bMuFitAny
            p=param(obj.muFitIndParam);
            mu(obj.bMuFit)=p(obj.muFitInd);
        end

        st =obj.stdFix;
        if obj.bStdFitAny
            p=param(obj.stdFitIndParam);
            st(obj.bStdFit)=p(obj.stdFitInd);
        end
        mustd=mu./st;

        cr =obj.crFix;
        if obj.bCrFitAny
            p=param(obj.crFitIndParam);
            cr(obj.bCrFit)=p(obj.crFitInd);
        end

        if nargin >= 3
            if size(rho,2) > 1
                obj.rho=rho(:,c);
            else
                obj.rho=rho;
            end
            if size(mu,2) > 1
                obj.mu =mu(:,c);
            else
                obj.mu =mu;
            end
            if size(cr,2) > 1
                obj.cr =cr(:,c);
            else
                obj.cr =cr;
            end
            if size(st,2) > 1
                obj.st =st(:,c);
            else
                obj.st =st;
            end
            if size(mustd,2) > 1
                obj.mustd =mustd(:,c);
            else
                obj.mustd =mustd;
            end
        else
            obj.rho=rho;
            obj.mu =mu;
            obj.cr =cr;
            obj.st =st;
            obj.mustd=mustd;
        end
    end
    function center_mu(obj)
        obj.muFixSet=obj.muFix;
        if ~obj.bCtrMu
            return
        end
        obj.muFixSet=repmat(obj.muFixSet,1,obj.ncmp);
        X=repmat((obj.cmpXunqC-obj.stdXunq)',numel(obj.muFix),1);
        ind=obj.muFixSet < 0;
        obj.muFixSet(ind)=X(ind);
        % TODO center bounds

    end
%- Bounds
    function get_bounds(obj)
        if ~isempty(obj.extType) && ~contains(obj.extType,'F')
            obj.get_ext_bounds();
            return
        else
            obj.lb=obj.LB;
            obj.ub=obj.UB;
        end

        % nans
        if any(isnan(obj.ub)) && ~any(isnan(obj.lb));
            [l,u]=obj.std_bound(obj.cmpXunqC);
            obj.lb(isnan(obj.lb))=l;
            obj.ub(isnan(obj.ub))=u;
        end
        %if ~isempty(obj.extType) && contains(obj.extType,'0')
        %    obj.lb(2)=0;
        %    obj.ub(2)=0;
        %end
    end

%- NLCon
    function get_nlcon(obj)
        obj.nlcon=[];

        if ~isempty(obj.extType)
            obj.get_ext_nlconst();
            return
            % HERE
        elseif obj.bObs
            obj.nlcon=@(p) obj.nlConstObs(p);
            %obj.nlcon=@(p) obj.nlcon_fun();
        end
    end
    function [c,ceq]=nlcon_fun(obj,param)
        c=[];
        ceq=[];

        % get rho
        rho=obj.rhoFix;
        if obj.bRhoFitAny
            p=param(obj.rhoFitIndParam);
            rho(obj.bRhoFit)=p(obj.rhoFitInd);
        end

        % positive definite
        c=DVFit.nlcon_posdef(rho);

        ceq=[];
    end
%- In0
    function get_In0(obj,bCombine)
        obj.In0=obj.F00;
        if ~bCombine
            return
        end

        obj.In0.CRL=repmat(obj.In0.CRL,1,1,obj.ncmp);
        obj.In0.CRU=repmat(obj.In0.CRU,1,1,obj.ncmp);
        if size(obj.In0.Mu,2)==1
            obj.In0.Mu=repmat(obj.In0.Mu,1,obj.ncmp);
        end
        obj.In0.Mu=permute(obj.In0.Mu,[3 1 2]);
    end
%- Out0
    function get_out0(obj)
        if obj.nSplit
            obj.out0=cell(1,obj.nSplit);
            for i = 1:obj.nSplit
                inds=(1:obj.nPass)+(i-1)*obj.nPass;
                obj.out0{i}=obj.get_out0_fun(inds);
            end
        else
            obj.out0=obj.get_out0_fun();
        end
    end
    function PP=get_out0_fun(obj,inds)
        % SELECT OUTPUT AND DATA
        if obj.bCombineCmp
            C=obj.nCmp;
        else
            C=1;
        end

        % split inds
        if nargin < 2
            inds=1:size(obj.rcmps{1},2);
        end

        PP=struct();

        PP.P=zeros(obj.nQuad,C);   % includes Pnn, Pnp...
        PP.N=zeros(obj.nQuad,C);   % includes Nnn, Nnp...
        PP.PA=zeros(obj.nAgree,C); % includes PD
        PP.PC=zeros(1,C);
        PP.nTrl=zeros(1,C);

        bAuto=size(obj.RCmpChs,3)==1;

        % HERE
        if obj.bBinary
            for c = 1:obj.ncmp
                for iq = 1:obj.nQuad

                    % PP.N is DATA THAT GOES INTO MODEL
                    if bAuto
                        r=obj.rcmps{c}(:,inds);
                        PP.N(iq,c)=sum(r==iq);
                    else
                        r1=obj.rcmps{c}(:,1);
                        r2=obj.rcmps{c}(:,2);
                        % XXX
                        r=r1==r2;
                        PP.N(iq,c)=sum(r==iq);
                    end

                    PP.nTrl(1,c)=numel(r);
                end
            end
        end
        PP.nTrlAll = sum(PP.nTrl(:));
    end
%- Param0
    function get_param0(obj)
        nError=1000;
        nBest=20;

        out=nan;
        c=0;
        best=inf;
        while true
            c=c+1;
            if c > nError
                error('Parameters do not appear to allow for valid evaluation')
            end
            p00=obj.param0_fun();
            out=obj.nllFun(p00);
            if isnan(out) || isinf(out)
                continue
            elseif out < best
                best=out;
                p0=p00;
            end
            if c > nBest && (~isnan(best) && ~isinf(best))
                break
            end
        end
        obj.param0=p0;
    end
    function param0=param0_fun(obj)
        if ~isempty(obj.extType) && ~contains(obj.extType,'F')
            param0=obj.get_ext_param0();
        else
            if obj.bCombineCmp
                val=1/8;
                C=obj.nCmp;
            else
                val=1/4;
                C=1;
            end
            param0  = DVFit.rand_interval_fun([obj.lb; obj.ub]).*val;
        end
    end
%- Obj
    function get_objective(obj)
        if ~isempty(obj.extType);
            obj.get_ext_objective();
            return
        end

        if obj.bBinary
            if obj.nPass==2
                f=@(p) obj.negLLFun_2_binary(p,obj.out0);
            else
                f=@(p) obj.negLLFun_n_binary(p,obj.out0);
            end
        else
            % TODO
        end
        obj.nllFun=f;
    end
%- Nll
    function [negLLAll,nLL,Out] = negLLFun_n_binary(obj,params,Out)
        obj.get_params(params);
        F=obj.insert_params_n(obj.In0,true);

        % Hadle bad rho
        flag=DVFit.checkRho(F.Rho);
        if ~isempty(flag)
            negLLAll=nan;
            nLL=nan;
            return
        end

        %x=tic;
        %Out.P=obj.autocdf(F,Out);
        %toc(x)

        Out.P=obj.manualcdf(F,Out);
        %diff(cat(3,Out.P,Out.P2),[],3)

        % OUT.P is the model
        % OUT.N is the data
        nLL=-sum(Out.N .* log(Out.P),1)/Out.nTrl;
        negLLAll = sum(nLL(:));

    end
%-cdf
    function P=autocdf(obj,F,Out)
        z=zeros(size(F.Mu(1,:,1)));
        P=Out.P;
        for c = 1:obj.ncmp
            % first has upperlimits
            P(1,c)=mvncdfcore(...
                            F.CRL(1,:,c),...
                            F.CRU(1,:,c),...
                            F.Rho, ...
                            1,...
                            obj.nPass,...
                            false,...
                            true,...
                            obj.mvnopts...
            );
            for i = 2:(obj.nQuad-1)
                P(i,c)=mvncdfcore(...
                                F.CRL(i,:,c),...
                                F.CRU(i,:,c),...
                                F.Rho, ...
                                1,...
                                obj.nPass,...
                                false,...
                                false,...
                                obj.mvnopts...
                );
            end

            %P(i,c)=mvncdf(...
            %                F.CRL(i,:,c),...
            %                F.CRU(i,:,c),...
            %                z,...
            %                F.Rho, ...
            %                'TolFun', 1e-12 ...
            %);

        end
        % last is remainder
        P(end,:)=1-sum(P,1);
    end
    function P=mvnl4(obj,F,Out);
        P=Out.P;
        for c = 1:obj.ncmp
            for i = 1:(obj.nQuad-1)
                P(i,c)=qvnl(...
                                F.CRU(i,:,c),...
                                F.Rho, ...
                                obj.mvnopts.TolFun ...
                );
            end
        end
        P(end,:)=1;
    end
    function [Q,P]=manualcdf(obj,F,Out)
        switch obj.nPass
        case 4
            P=obj.mvnl4(F,Out);
        otherwise
            TODO
        end
        Q=Out.P;

        % first
        Q(1,:)=P(1,:);

        for c = 1:obj.ncmp
            Ind=F.CRU(:,:,c)~=inf;
            M=sum(Ind,2); % number of crits
            for ki =1:obj.nPass-1
                k=obj.nPass-ki;
                cind=find(M==k);

                for l = 1:length(cind)
                    i=cind(l);
                    ind=Ind(i,:);

                    % ki = num infs
                    % k  = num crits
                    % skips ki={0,last}
                    if ki==1
                        Q(i,c)=P(i,c)-Q(1,c);
                    elseif ki==2
                        % matching overlapping crit
                        % matches 2 for nPass==4
                        bN=any(bsxfun(@times,~Ind,~ind),2); % matching crit pattern
                        n1=find(bN & M ==(k+1));

                        Q(i,c)=P(i,c)-sum(P(n1,c))+Q(1,c);

                    elseif ki==3
                        %bN=any(bsxfun(@times,~Ind,~ind),2);
                        %n1=find(bN & M ==(k+1));
                        %n=P(i,c)-sum(P(n1,c))+sum(Q(m2,c))+2*Q(1,c);

                        % matching overlapping inf
                        % matches 3 for nPass=3
                        bM=any(bsxfun(@times,Ind,ind),2);  % matching inf pattern
                        m1=find(bM & M ==(k+1));
                        m2=find(bM & M ==(k+2)); % same as n2 (always ???)

                        m=P(i,c)-sum(P(m1,c))+sum(Q(m2,c))+2*Q(1,c);

                        Q(i,c)=P(i,c)-sum(P(m1,c))+sum(Q(m2,c))+2*Q(1,c);
                    end

                end
            end
        end
        % last
        Q(end,:)=1-sum(Q(1:end-1,:),1);
    end
%- insert
    function F=insert_params_n(obj,F,bInit)
        % Rho
        if bInit || obj.bRhoFitAny
            F.Rho(F.indRhoU)=obj.rho;
            if F.bRhoD
                L=F.Rho';
                L(F.teye)=0;
                F.Rho=F.Rho + L;
            else
                F.Rho=F.Rho + F.Rho' + F.eye;
            end
        end

        if bInit || obj.bMuCrFitAny
            % STD ALREADY HANDLED

            % Mu
            if bInit || obj.bMuFitAny
                F.Mu(:)=obj.mustd(:);
            end

            % Cr
            if bInit || obj.bCrFitAny
                for i = 1:obj.nPass
                    F.CRL(F.CRL==i)=obj.cr(i);
                    F.CRU(F.CRU==i)=obj.cr(i);

                    F.CRB(F.CRB==i)=obj.cr(i);
                end
            end

            % combine mu and cr
            F.CRL=bsxfun(@minus,F.CRL,F.Mu);
            F.CRU=bsxfun(@minus,F.CRU,F.Mu);

            F.CRB=bsxfun(@minus,F.CRB,F.Mu);
        end

    end
%- Agree
    function Out=getAgree(obj,Out)
        % PA
        for ia = 1:obj.nAgree
            Out.PA(ia,:)=sum( Out.P(obj.agreeInds==ia,:) ,1);
        end
    end
    function Out=getPCmp(obj,Out)
        % PC
        Out.PC=sum(Out.P./(obj.nAgree./obj.agreeInds),1);
    end

end
methods(Static)
    function flag=checkRho(rho)
        ei=eig(rho);
        if any(ei < 0)
            flag=-1;
        elseif any(ei < 10^-10)
            flag=0;
        elseif any(abs(rho) > 1)
            flag=2;
        else
            flag=[];
        end
    end
    function r=rand_interval_fun(range)
        if size(range,1)==1 && size(range,2)==2
            range = transpose(range);
        end
        m=size(range,2);

        rng('shuffle');
        rnd=rand(1,m);

        d=diff(range,[],1);
        r=(rnd.*d)-d/2;
    end
    function [L,U]=std_bound(cmpX)
        % m0 s0 b0
        %L = [2.0.*min(cmpX-mean(cmpX))+mean(cmpX) 0.02.*(max(cmpX)-min(cmpX)) 0.35];
        %U = [2.0.*max(cmpX-mean(cmpX))+mean(cmpX) 2.00.*(max(cmpX)-min(cmpX)) 3.00];
        L = [0.02.*(max(cmpX)-min(cmpX))];
        U = [2.00.*(max(cmpX)-min(cmpX))];
    end
    function c=nlcon_posdef(rho)
        % keep rho positive definite
        % https://www.mathworks.com/matlabcentral/answers/632024-imposing-positive-semi-definiteness-and-symmetry-in-fmincon
        %
        c=-eig(diag(1-exp(-10*rho)));
        % positive-definintes
        %c=-DVFit.texp(rho);
        %c=-DVFit.texp(rho);
    end
end
end

