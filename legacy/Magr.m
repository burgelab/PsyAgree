classdef Magr < handle & Magr_plot
properties
    PEmp %
    PFit  % Fit
    PSim  % SIM
    PPoly % POLY
    PBino % BINOMIAL
    PBoot % BOOT & old Emp
end
properties(Hidden)
    mvnopts
    G
end
methods
    function obj=Magr(Data)
        %obj.mvnopts=RMvnCdf.getMVNOpts();
        %obj.G=RMvnCdf.getGStruct();
    end
    function obj=get_magr(obj,nSim,bProg)
        if nargin < 2; nSim=[]; end
        if nargin < 3; bProg=false; end


        % NO SIMULATION
        obj.get_agree_emp();
        obj.get_agree_boot();

        obj.get_agree_fit(nSim); % CURVE
        obj.get_agree_binom(nSim,bProg);
        obj.get_agree_sim(nSim,bProg);

        % FITTED PARAB FROM BOOT/EMP
        % obj.get_agree_parab();
    end
    function get_agree_emp(obj)
        C=numel(obj.cmpXunq);
        obj.PEmp.PA=zeros(C,1);
        obj.PEmp.PC=zeros(C,1);
        cmpXunq=unique(obj.CMPX);
        for c = 1:C
            IND=all(obj.CMPX==cmpXunq(c),2);
            R=obj.RCMPCHS(IND,:);
            obj.PEmp.PC(c)=sum(R(:))./numel(R);
            obj.PEmp.PA(c)=sum(sum(diff(R,[],2),2)==0,1)/size(R,1);
        end
    end
    function get_agree_boot(obj,nBoot)
        C=numel(obj.cmpXunq);
        if nargin < 2 || isempty(nBoot)
            nBoot=obj.DVFITTER.nBoot;
            if nBoot == 1
                nBoot=100;
            end
        end

        PBoot.PC=zeros(C,nBoot);
        PBoot.PA=zeros(C,nBoot);
        for c = 1:C
            IND=all(obj.CMPX==obj.cmpXunq(c),2);
            R=obj.RCMPCHS(IND,:);
            for b = 1:nBoot
                N=size(R,1);
                ii=randi(N,N,1);
                r=R(ii,:);

                PBoot.PC(c,b)=sum(r(:))./numel(r);
                PBoot.PA(c,b)=sum(sum(diff(r,[],2),2)==0,1)/size(r,1);
            end
        end
        PAm=mean(PBoot.PA,2);
        PCm=mean(PBoot.PC,2);
        obj.PBoot.PA=PBoot.PA;
        obj.PBoot.PC=PBoot.PC;
        obj.PBoot.PAm=PAm;
        obj.PBoot.PCm=PCm;
    end
%- AGREE CURVE
    function obj=get_agree_fit(obj,nSim,bProg)
        if numel(unique(obj.FIT.RHO))==1
            rho=obj.FIT.RHO(1);
        else
            rho=mean(obj.FIT.RHO);
        end
        % CURVE
        [obj.PFit.PA,obj.PFit.PC]=obj.getAgreeCurve(rho,[],[],obj.FIT.CR1,obj.FIT.CR2,obj.G,obj.mvnopts);
    end
    function obj=get_agree_sim(obj,nSim,bProg)
        if nargin < 2; nSim=[];  end
        if nargin < 3; bProg=[]; end

        if numel(unique(obj.FIT.RHO))==1
            rho=obj.FIT.RHO(1);
        else
            rho=mean(obj.FIT.RHO);
        end
        mu=-2:0.2:2;
        str='Fitting simulations';
        [obj.PSim.PA,obj.PSim.PC,obj.PSim.PA68L,obj.PSim.PA68U,obj.PSim.PA95L,obj.PSim.PA95U]=obj.agree_curve(rho,mu,nSim,bProg,str);
    end
    function obj=get_agree_binom(obj,nSim,bProg)
        if nargin < 2; nSim=[];  end
        if nargin < 3; bProg=[]; end
        rho=0;
        mu=-2:0.2:2;
        str='Fitting binomial model';
        [obj.PBino.PA,obj.PBino.PC,obj.PBino.PA68L,obj.PBino.PA68U,obj.PBino.PA95L,obj.PBino.PA95U]=obj.agree_curve(rho,mu,nSim,bProg,str);
        %obj.PBino.PA
    end
    function obj=get_agree_parab(obj)
        obj.get_magr();

        if obj.nBoot<=1
            PC=obj.PEmp.PC;
            PA=obj.PEmp.PA;
        else
            PC=obj.PBoot.PC;
            PA=obj.PBoot.PA;
        end
        order=10;

        CI = [0.5.*(1-obj.CI/100) 1-0.5.*(1-obj.CI/100)];
        xfix=[-.5 .5];
        yfix=[1 1];


        nBoot=size(PC,2);
        p=cell(nBoot,1);
        o=order;
        while true
            for i = 1:nBoot
                p{i}=polyfix(PC(:,i)-.5,PA(:,i),2,xfix,yfix);
            end
            sz=cellfun(@numel,p);
            if all(diff(sz)==0)
                o=sz(1);
                p=vertcat(p{:});
                break
            else
                b=unique(sz);
                [count]=hist(sz,b);
                o=b(b==max(count));
            end
        end
        Uyi=quantile(p(:,end),CI(1));
        Lyi=quantile(p(:,end),CI(2));
        p=zeros(nBoot,o);
        for i = 1:nBoot
            pU(i,:)=polyfix(PC(:,i)-.5,PA(:,i),o-1,[xfix 0],[yfix Uyi]);
            pL(i,:)=polyfix(PC(:,i)-.5,PA(:,i),o-1,[xfix 0],[yfix Lyi]);
        end

        U=mean(pU,1);
        L=mean(pL,1);
        m=mean([U;L],1);

        obj.PPoly.PC=m;
        obj.PPoly.PCU=U;
        obj.PPoly.PCL=L;
    end
    function [PA,PC,PA68l,PA68u,PA95l,PA95u]=agree_curve(obj,rho,mu,nSim,bProg,str)
        if nargin < 2 || isempty(rho)
            rho=obj.FIT.RHO;
        end
        if nargin < 3 || isempty(mu)
            mu=-2:.2:2;
        end
        if nargin < 4 || isempty(nSim)
            nSim=obj.nParabSim;
        end
        if nargin < 5 || isempty(bProg)
            bProg=true;
        end
        if isempty(obj.nTrlParabSim) || obj.nTrlParabSim==0
            nTrlParabSim=size(obj.CMPX,1);
        else
            nTrlParabSim=obj.nTrlParabSim;
        end

        COV=[1 rho; rho 1];
        mu=repmat(transpose((mu)),[1 2]);
        CI68 = [0.5.*(1-68.2/100) 1-0.5.*(1-68.2/100)];
        CI95 = [0.5.*(1-95.4/100) 1-0.5.*(1-95.4/100)];

        flag=0;
        N=2; % XXX total n passes
        magr=2; % XXX how many of n passes we want to agree
        if isequal(size(COV),[N N])
            cov=COV;
        elseif size(COV,1)==size(mu,1) && size(COV,2)==1
            flag=1;
        end

        if bProg
            p=Pr(length(mu),length(mu)/100,str);
        end
        for i = 1:length(mu) % FOR EACH MEAN (cmp)
            if flag
                cov=COV(i)*eye(N);
                cov(cov==0)=1;
                cov=flip(cov,2);
            end
            if bProg
                p.u();
            end


            % COMPUTE AGREEMENT AND CMP CHOSEN MONTE CARLO SIMULATIONS
            [PAtmp, PCtmp] = obj.get_agree_curve_mc(N,magr,mu(i,:),cov,nSim,nTrlParabSim);
            % STORE ACTUAL P(A) AND P(C)
            PC(i,:) = mean(PCtmp);
            PA(i,:) = mean(PAtmp);
            % STORE LOWER BOUNDS ON P(A)
            PA68l(i,:) = quantile(PAtmp,CI68(1));
            PA68u(i,:) = quantile(PAtmp,CI68(2));
            PA95l(i,:) = quantile(PAtmp,CI95(1));
            PA95u(i,:) = quantile(PAtmp,CI95(2));
        end
        if bProg
            p.c();
        end
    end
    function [PA,PC]=get_agree_curve_mc(obj,N,M,mu,cov,numSim,numTrlPerSim)
    % simulates 2AFC trials for which the observer's decision variable is
    % Gaussian Set.distributed, then computes proportion agreement and proportion
    % cmp chosen
        if sum(cov(logical(1-eye(N))))==0 % IF THERE IS NO CORRELATION BETWEEN DECISION VARIABLES, CAN COMPUTE FAST
            % GENERATE DECISION VARIABLES
            Z = normrnd(mu(1),cov(1,1),[numTrlPerSim N numSim]);
            % COMPUTE PROPORTION CMP CHOSEN
            PC = transpose(sum(reshape(Z,[size(Z,1)*size(Z,2) size(Z,3)])>0))./(numTrlPerSim*N);
            % COMPUTE AGREEMENT
            PA = squeeze(sum(sum(Z>=0,2)==M | sum(Z<0,2)==M)./numTrlPerSim);
        else % IF THERE IS CORRELATION BETWEEN DECISION VARIABLES, LOOP (SLOWER)
            for i = 1:numSim % FOR EACH SIMULATION
                % GENERATE DECISION VARIABLES
                Z = mvnrnd(mu,cov,numTrlPerSim);
                % COMPUTE PC
                PC(i,:) = sum(Z(:)>0)./numel(Z);
                % COMPUTE M-AGREEMENT
                bZmoreThan0 = Z>=0;
                bZlessThan0 = Z<0;
                PA(i,:) = sum(sum(bZmoreThan0,2)==M | sum(bZlessThan0,2)==M)./numTrlPerSim;
            end
        end
    end
    function [x,y,X,Y]=get_agree_interp(obj,PC,PA,PAL,PAU)
        n=numel(obj.RCMPCHS);
        PCbinRes = 1./n;
        PCbinInterp   = transpose([PCbinRes:PCbinRes:1-PCbinRes]);
        PAbinInterp   = pchip(PC,PA, PCbinInterp);
        PAl           = pchip(PC,PAL,PCbinInterp);
        PAu           = pchip(PC,PAU,PCbinInterp);
        PCi           = PCbinInterp;

        %constraint = [1,0,1,1];  % Ignore possible x^2 term
        %polyFunc = @(p) polyval(p.*constraint,x);
        %objectiveFunc = @(p) (y - polyFunc(p)).^2;
        %p0 = [1, 0, 2, 3];  % It pays to have a realistic initial guess
        %p = fminsearch( objectiveFunc, p0 );

        x=PCbinInterp;
        y=PAbinInterp;

        X=[PCi; flipud(PCi)];
        Y=[PAl; PAu];
    end
end
methods(Static=true)
    function [pA,pC]=getAgreeCurve(rho,mu1,mu2,cr1,cr2,G,mvnopts)
        if nargin < 2 || isempty(mu1)
            mu1=[-4:.1:4];
        end
        if nargin < 3 || isempty(mu2)
            mu2=mu1;
        end
        if nargin < 4 || isempty(cr1)
            cr1=0;
        end
        if nargin < 5 || isempty(cr2)
            cr2=cr1;
        end
        if nargin < 6 || isempty(G)
            G=RMvnCdf.getGStruct();
        end
        if nargin < 7 || isempty(mvnopts)
            mvnopts=RMvnCdf.getMVNOpts();
        end


        [Rho]=Magr.parseRho(rho);
        [CRL,CRU]=Magr.getCR(cr1(1),cr2(2));
        nDim=size(CRL,2);
        cls=class(Rho);

        pA=zeros(length(mu1),1);
        pC=zeros(length(mu1),1);
        for i=1:numel(mu1)
           MU=Magr.parseMu(mu1(i),mu2(i));
           [CRL0,CRU0]=Magr.getCR0(CRL,CRU,MU,[1 1],nDim);
           p1=prod(RMvnCdf.Phi(CRU0),2); % XXX SLOW 1
           if numel(Rho)>1
               Rho=Rho(2);
           end
           %[Ppp,Pnn,Ppn,Pnp,pC(i),pA(i)]=Magr.fitd_r_pat(CRL0,CRU0,Rho,cls,G,p1,mvnopts);
           [Ppp,Pnn,Ppn,Pnp,pC,pA]=Magr.fitd_r_pat(CRL0,CRU0,Rho,cls,G,p1,mvnopts);


           %[Ppp,Pnn,Ppn,Pnp,pC(i),pA(i)]=...
           %    Magr.getFittedResponsePatternOld(rho,mu1(i),mu2(i),cr1(1),cr2(1),mvnopts);
        end
    end
    function [varE,varI,varT]=getNoiseSourcesInd(cmp,std,rho,dp)
        if numel(unique(cmp))==1 && numel(unique(std))==1
            cmp=cmp(1);
            std=std(1);
        end
        varT=((cmp-std)/dp).^2;
        varE=rho*varT;
        varI=varT-varE;
    end
    function [Pnn,Pnp,Ppn,Ppp,pC,pA,pD]=getFittedResponsePattern(rho,mu1,mu2,cr1,cr2,G,mvnopts)
        if nargin < 6 || isempty(G)
            G=RMvnCdf.getGStruct();
        end
        if nargin < 7 || isempty(mvnopts)
            mvnopts=RMvnCdf.getMVNOpts();
        end
        % RMvnCdf is modified version of mvncdf

        MU=Magr.parseMu(mu1,mu2);
        [SIG,s,Rho]=Magr.parseRho(rho);

        [CRL,CRU]=Magr.getCR(cr1,cr2);
        [CRL0,CRU0]=Magr.getCR0(CRL,CRU,MU,s);
        cls=class(Rho);

        [Pnn,Pnp,Ppn,Ppp,pC,pA,pD]=...
            Magr.fitd_r_pat(CRL0,CRU0,Rho(2),cls,G,mvnopts);
    end
    function [Pnn,Pnp,Ppn,Ppp,pC,pA,pD]=getFittedResponsePatternOld(rho,mu1,mu2,cr1,cr2,mvnopts)
        Pnn = mvncdf([cr1 cr2],[mu1 mu2],[1 rho; rho 1]);        % NEG/NEG
        Pnp = mvncdf([cr1 Inf],[mu1 mu2],[1 rho; rho 1]) - Pnn;  % NEG/POS
        Ppn = mvncdf([Inf cr2],[mu1 mu2],[1 rho; rho 1]) - Pnn;  % POS/NEG
        Ppp = 1 - Pnn - Ppn - Pnp;                               % POS/POS
        pC=Ppp + Ppn/2 + Pnp/2;
        pA=Ppp + Pnn;
        pD=Ppn + Pnp;
    end
end
methods(Static,Hidden)
    function MU=parseMu(mu1,mu2)
        % DO ONCE PER UNIQUE MU
        MU =[mu1 mu2];
    end
    function [Rho,s]=parseRho(rho,bSkip)
        if nargin < 2 || isempty(bSkip)
            bSkip=false;
        end
        % DO ONCE PER UNIQUE SIGMA
        Rho=[1 rho; rho 1];
        s=[1 1];
        if ~bSkip
            [~,err] = cholcov(Rho,0); % XXX SLOW
            if err ~= 0
                error(message('stats:mvncdf:BadMatrixSigma'));
            end
        end
        %s = sqrt(diag(Sig))';
        %Rho = Sig ./ (s'*s);
    end
    function [CRL,CRU]=getCR(cr1,cr2,cr0)
        % DO ONCE PER UNIQUE cr1 & cr2
        if nargin < 3 || isempty(cr0)
            cr0=Inf(3,2);
        end
        CRL=cr0;
        CRU=cr0;
        CRU(1,:)=[cr1 cr2]; % XXX SLOW 3
        CRU(2,1)=cr1;
        CRU(3,2)=cr2;
    end
    function [CRL0,CRU0]=getCR0(CRL,CRU,MU,s,nDim)
        if nargin < 5 || isempty(nDim)
            nDim=size(CRL,2);
        end
        % DO ONCE PER UNIQUE cr1, cr2, mu1, mu2, rho
        %if nargin < 6
        %    c0=cell(3,1);
        %end
        %CRL0=c0;
        %CRU0=c0;
        if nDim > 2
            CRL0=(CRL-MU)./s;
            %for i = 1:3
                %CRL0{i}=cr0_fun(CRL{i},MU,s);
            %end
        else
            CRL0=CRL;
        end
        CRU0=(CRU-MU)./s;
        %for i = 1:3
            %CRU0{i}=(CRU{i}-MU)./s;
        %end

        %CRU0{1}=cr0_fun(CRU{1},MU,s);
        %CRU0{2}=cr0_fun(CRU{2},MU,s);
        %CRU0{3}=cr0_fun(CRU{3},MU,s);

        function [CR0]=cr0_fun(CR,MU,s)
            % XXX make sure this is needed
            CR0 = bsxfun(@rdivide, bsxfun(@minus,CR,MU),s);
        end
    end
    function [Pnn,Pnp,Ppn,Ppp,pC,pA,pD]= fitd_r_pat_old(CRL0,CRU0,MU,SIG,s,Rho,cls,G,mvnopts)

        Pnn = RMvnCdf.get(CRL0{1},CRU0{1},MU,SIG,s,Rho,cls,G,mvnopts);        % NEG/NEG
        Pnp = RMvnCdf.get(CRL0{2},CRU0{2},MU,SIG,s,Rho,cls,G,mvnopts) - Pnn;  % NEG/POS
        Ppn = RMvnCdf.get(CRL0{3},CRU0{3},MU,SIG,s,Rho,cls,G,mvnopts) - Pnn;  % POS/NEG

        %Pnn = RMvnCdf.get(CRL0{1},CRU0{1},MU,SIG,s,Rho,cls,G,mvnopts);        % NEG/NEG
        %Pnp = RMvnCdf.get(CRL0{2},CRU0{2},MU,SIG,s,Rho,cls,G,mvnopts) - Pnn;  % NEG/POS
        %Ppn = RMvnCdf.get(CRL0{3},CRU0{3},MU,SIG,s,Rho,cls,G,mvnopts) - Pnn;  % POS/NEG

        Ppp = 1 - Pnn - Ppn - Pnp;                      % POS/POS
        pC=Ppp + Ppn/2 + Pnp/2;
        pA=Ppp + Pnn;
        pD=Ppn + Pnp;
    end
    function [pP,pA,pC]= fitd_r_pat(CRL0,CRU0,Rho,mvnopts,mvnconsts,agreeInds,phi)
        % F.p1 = prod(0.5 * erfc(-F.CRU0 / sqrt(2)),2); % PHI

        nA=max(agreeInds);

        % init
        pP=zeros(mvnconsts.nQuad,1,mvnconsts.cls);
        pA=zeros(nA,1,mvnconsts.cls);
        pC=zeros(1,1,mvnconsts.cls);

        if mvnconsts.nDim==2
            % Pnn Pnp Ppn Ppp
            mvnconsts=RMvnCdf.critG(Rho,mvnconsts);
            for i = 1:3
                pP(i) = RMvnCdf.bvncdf(...
                    CRU0(i,:),...
                    Rho,...
                    mvnconsts.cls,...
                    mvnconsts.w,...
                    mvnconsts.y,...
                    mvnconsts.p2,...
                    phi(i),...
                    true ...
                );
            end
            pP(end) = 1 - sum(P);

            %pA(1)=pP(4) + pP(1); % agree
            %pA(2)=pP(2) + pP(3); % disagree

            %pC=pP(4) + pP(2)/2 + pP(3)/2;
        else
            for i = 1:mvnconsts.nQuad
                % XXX faster general method?
                pP(i)=RMvnCdf.get(CRL0(i,:),CRU0(i,:),Rho,mvnconsts);
            end

        end
        for ia = 1:nA
            pA(i)=sum(pP(agreeInds==ia));
        end
        pC=pP./(nA./obj.agreeInds);
    end
    function [patterns,agreeInds]=getRespPats(R,nPass)
    % get response patterns
        rs=unique(R(:));
        pats=repmat({rs},1,nPass);
        pats=Set.distribute(pats{:});
        [~,ind]=sort(sum(pats,2));
        patterns=pats(ind,:);
        [~,~,agreeInds]=unique(sum(patterns,2));
    end
end
end
