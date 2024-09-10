classdef DVExt < handle
properties(Hidden)
    Fake

    bExt=0
    bFake=0
    bZero=0

end
methods
    function F=get_params_obs(obj,p)
        F=struct();

        obj.get_params(p);
        F.stdT1=obj.stdFix(1);
        F.stdT2=obj.stdFix(4);

        F.rho1=obj.rho(1);
        F.rho12=obj.rho(2);
        F.rho2=obj.rho(end);

        F.varS=F.rho12 * F.stdT1 * F.stdT2;
        F.varE1=F.rho1 * F.stdT1^2;
        F.varE2=F.rho2 * F.stdT2^2;
        F.varI1=F.stdT1^2-F.varE1;
        F.varI2=F.stdT2^2-F.varE2;
        %F.varP1=F.varE1-F.varS;
        %F.varP2=F.varE2-F.varS;

        F.Rho=p;
    end
    function [p1,p2]=get_params_split(obj,p,c)
        p1=[p(2:3) p(1)];
        p2=[p(4:5) p(1)];

        if nargin >= 3
            r=cell(1,2);
            for i = 1:obj.nSplit
                obj.get_params(p1,c);
                r{i}=obj.rho;
            end
        end
    end
    function F=get_params_ext(obj,p)
        % p(1) stdB
        % p(2) covLB
        % p(3) stdL
        F=struct();

        F.stdB=obj.stdFix(1);
        F.stdL=obj.stdFix(4);

        if length(p)==2
            F.varB=p(1);
            F.covLB=0;
            F.varL=p(2);
        else
            F.varB=p(1);
            F.covLB=p(2);
            F.varL=p(3);
        end

        [F.rB,F.rBL,F.rL,F.corrLB]=Ext.vars2rhos(F.varB,F.covLB,F.varL,F.stdB,F.stdL);


        F.Rho=[F.rB; F.rBL; F.rL];
    end
    function F=get_params_ext_fake(obj,p)
        obj.get_params(p);

        F=struct();
        F.stdB=obj.stdFix(1);
        F.stdL=obj.stdFix(4);
        [F.varB,F.varL,F.covLB,F.corrLB,F.varIB,F.varIL,F.varEB]=Ext.rhos2vars(obj.rho(1),obj.rho(2),obj.rho(3),F.stdB,F.stdL);
    end
%- Bounds
    function get_ext_bounds(obj)
        switch obj.extType
        case {'E','E0'}
            s=obj.stdFix;
            varB=s(1).^2;
            varL=s(4).^2;
            covLB=(sqrt(varB)*sqrt(varL));

            L=zeros(3,1);
            U=zeros(3,1);
            L(1:3)=0.000001*3600;

            U(3)=varL;
            U(2)=covLB-L(3);
            L(2)=-U(2);
            U(1)=varB -L(3)-2*L(2);

            nl=@(p) obj.nlContExt(p);
            if contains(obj.extType,'0')
                L(2)=[];
                U(2)=[];
            end
        case 'EL'
            TODO
        end
        obj.lb=L;
        obj.ub=U;
    end
%- Param0
    function param0=get_ext_param0(obj)
        switch obj.extType
        case {'E','E0'}
            rng=(obj.ub-obj.lb);
            param0 = rand(numel(obj.ub),1).*rng+obj.lb;
        case 'EL'
        end
    end
%- Objective
    function get_ext_objective(obj)
        obj.bFake=contains(obj.extType,'F');
        obj.bExt=~isempty(obj.extType) & ~obj.bFake;
        obj.bZero=contains(obj.extType,'0');

        switch obj.extType
        case {'E','E0'}
            f=@(p) obj.negLLFun_ext(p,obj.out0);
        case 'EL'
            f=@(p) obj.negLLFun_extL(p,obj.out0);
        case {'F','F0'}
            f=@(p) obj.negLLFun_extFake(p,obj.out0);
        end
        obj.nllFun=f;
    end
%- NegLL
    function [negLLAll,nLL]=negLLFun_ext_split(obj,p,PP1,PP2)
        [p1,p2]=obj.get_params_split(p);

        [negLLAll1,nLL1,PP1]=obj.negLLFun_extFake(p1,PP1);
        [negLLAll2,nLL2,PP2]=obj.negLLFun_extFake(p2,PP2);

        negLLAll=negLLAll1+negLLAll2;
        negLL=negLL1+negLL2;

    end
    function [negLLAll,nLL]=negLLFun_extFake_split(obj,p,PP1,PP2)
        [p1,p2]=obj.get_params_split(p);

        [negLLAll1,nLL1,PP1]=obj.negLLFun_extFake(p1,PP1);
        [negLLAll2,nLL2,PP2]=obj.negLLFun_extFake(p2,PP2);

        negLLAll=negLLAll1+negLLAll2;
        negLL=negLL1+negLL2;

    end
    function [negLLAll,nLL,PP]=negLLFun_ext(obj,p,PP)
        F=obj.get_params_ext(p);
        if any(abs(F.Rho)> 1)
            negLLAll=nan;
            negLL=nan;
            return
        end

        [negLLAll,nLL,PP]=obj.negLLFun_n_binary(F.Rho,PP);
    end
    function [negLLAll,nLL,Out]=negLLFun_extFake(obj,p,Out)
        Fk=obj.get_params_ext_fake(p);
        F=obj.insert_params_n(obj.In0,true);

        % Hadle bad rho
        flag=DVFit.checkRho(F.Rho);
        if ~isempty(flag)
            negLLAll=nan;
            negLL=nan;
            return
        end

        Out.P=obj.autocdf(F,Out);

        nLL=-sum(Out.N .* log(Out.P),1)/Out.nTrl;
        negLLAll = sum(nLL(:));
    end
    function [negLLAll,nLL,PP]=negLLFun_ext2(obj,params,PP)
        r=zeros(3,1);

        [negLLAll,nLL,PP]=obj.negLLFun_n_binary(r,PP);
    end
    function negLLFun_ext_linear(obj,params,PP)
        T1=B1_1*log(x)+B0_1;
        T2=B1_2*log(x)+B0_2;
    end
%- Constr
    function get_ext_nlconst(obj)
        switch obj.extType
        case {'F','F0'}
            obj.nlcon=[];
            obj.nlcon=@(p) obj.nlConstExtFake(p);
        otherwise
            obj.nlcon=@(p) obj.nlConstExt(p);
        end
    end
    function [C,Ceq]=nlConstExt(obj,p)
        F=obj.get_params_ext(p);

        C=[];
        Ceq=[];
        C=DVFit.nlcon_rho(C,F);
        %C=DVFit.nlcon_varB(C,F);
        %C=DVFit.nlcon_covLB(C,F);
        if contains(obj.extType,'0')
            Ceq=DVFit.nlcon_covLB0(Ceq,F);
        else
            C=DVFit.nlcon_corrLB(C,F);
        end

    end
    function [C,Ceq]=nlConstObs(obj,p)
        C=[];
        Ceq=[];
        F=obj.get_params_obs(p);
        C=DVFit.nlcon_obs(C,F);
    end
    function [C,Ceq]=nlConstExtFake(obj,p)
        Fk=obj.get_params_ext_fake(p);

        C=[];
        Ceq=[];
        %C=DVFit.nlcon_rho(C,Fk);
        %C=DVFit.nlcon_varB(C,Fk);
        %C=DVFit.nlcon_covLB(C,Fk);
        if contains(obj.extType,'0')
            Ceq=DVFit.nlcon_covLB0(Ceq,Fk);
        end

    end
end
methods(Static)
    function C=nlcon_rho(C,F)
        lt=1;
        gt=0;
        C=[C; -F.Rho(:)+gt; F.Rho(:)-lt];
        %C=[C; -F.Rho(:)+gt; F.Rho(:)-lt];
        %rho=[-F.Rho(:)-1; -F.corrLB(:)-1; F.Rho(:); F.corrLB(:)-1];
        %rho=[];
        %C=[C; -rho-1; rho-1];
        %C=[C; -rho-1; rho];
    end
    function C=nlcon_corrLB(C,F)
        lt=1;
        gt=-1;
        C=[C; -F.corrLB(:)+gt; F.corrLB(:)-lt];
    end
    function C=nlcon_varB(C,F)
        C=[C; [(F.varL + F.varB + F.covLB) - F.stdB.^2]];
    end
    function C=nlcon_covLB(C,F)
        C=[C; [F.covLB-abs(F.stdL .* F.stdB)]];
    end
    %
    function Ceq=nlcon_covLB0(Ceq,F)
        Ceq=[Ceq; F.covLB];
    end
    function C=nlcon_obs(C,F)

        gt=0;
        C=[
            %-F.varP1+gt;
            %-F.varP2+gt;
            -F.varE1+gt;
            -F.varE2+gt;
            -F.varI1+gt;
            -F.varI2+gt;
             %F.varS-F.varE1;
             %F.varS-F.varE2;
             F.varS-sqrt(F.varE2.*F.varE1);
        ];

    end

end
end
