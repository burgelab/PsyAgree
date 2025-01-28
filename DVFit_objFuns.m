classdef DVFit_objFuns < handle
properties
end
methods(Static)
end
methods
%- main
    function [negLLAll,nLL,Out] = negLLFun_2_binary(obj,params,Out)
        obj.get_params(params);
        %F=obj.insert_params_2(obj.F0,true);
        F=obj.insert_params_n(obj.In0,true);
        rho=F.Rho(2);

        % consts
        C=RMvnCdf.critG(rho,obj.mvnconsts);
        phi = prod(0.5 * erfc(-F.CRB./sqrt(2)), 2);

        for c = 1:obj.ncmp
            Out.P(1,c) = RMvnCdf.bvncdf(F.CRB(1,:,c),rho,C.cls, C.w, C.y, C.p2, phi(1,:,c),true);
            Out.P(2,c) = RMvnCdf.bvncdf(F.CRB(2,:,c),rho,C.cls, C.w, C.y, C.p2, phi(2,:,c),true)-Out.P(1,c);
            Out.P(3,c) = RMvnCdf.bvncdf(F.CRB(3,:,c),rho,C.cls, C.w, C.y, C.p2, phi(3,:,c),true)-Out.P(1,c);
        end
        Out.P(4,:) = 1-sum(Out.P,1);

        % PA
        for ia = 1:obj.nAgree
            Out.PA(ia,:)=sum( Out.P(obj.agreeInds==ia,:) ,1);
        end

        % PC
        Out.PC=sum(Out.P./(obj.nAgree./obj.agreeInds),1);


        % NEGATIVE LOG-LIKELIHOOD PER DATA POINT
        nLL=-sum(Out.N .* log(Out.P),1); % over pattern type
        negLLAll = sum(nLL(~isnan(nLL(:)))); %over comparisons

        if obj.bStdFitAny
            nLLS=obj.negLLStd();
            negLLAll=negLLAll+nLLS;
        end

    end
%- CDF
    function P=fullcdf(obj,F,Out);
        P=Out.P;
        for c = 1:obj.ncmp
            for i = 1:(obj.nQuad)
                P(i,c)=qvncdf(...
                                F.CRU(i,:,c),...
                                F.Rho, ...
                                'TolFun', 1e-12 ...
                );
            end
            P(end,c)=1;
        end
    end
    function [Q,P]=manualcdf(obj,F,Out)
        P=obj.fullcdf(F,Out);
        Q=Out.P;

        for c = 1:obj.ncmp
            % first
            Q(1,c)=P(1,c);

            Ind=F.CRU(:,:,c)~=inf;
            M=sum(Ind,2);
            for ki =1:obj.nPass-1
                k=obj.nPass-ki;
                cind=find(M==k);

                for l = 1:length(cind)
                    i=cind(l);
                    ind=Ind(i,:);

                    if k==3
                        Q(i,c)=P(i,c)-Q(1,c);
                        continue;
                    end

                    %m1
                    bM=any(bsxfun(@times,~Ind,~ind),2);
                    m1=find(bM & M ==(k+1));
                    if k==2
                        Q(i,c)=P(i,c)-sum(P(m1,c))+Q(1,c);
                        continue
                    end

                    bM=any(bsxfun(@times,Ind,ind),2);
                    m1=find(bM & M ==(k+1));
                    bM=any(bsxfun(@times,Ind,ind),2);
                    m2=find(bM & M ==(k+2));
                    %if all(ind == [ 1 0 0 0 ])
                    %    ind
                    %    Ind(m1,:)
                    %    Ind(m2,:)
                    %end


                    %if i==15
                    %    Ind(m1,:)
                    %    Ind(m2,:)
                    %    ind
                    %elseif k==1
                        Q(i,c)=P(i,c)-sum(P(m1,c))+sum(Q(m2,c))+2*Q(1,c);
                    %end
                end
            end

            % last
            Q(end,c)=1-sum(Q(:,c),1);
        end
    end
%- STd
    function negLL=negLLStd(obj)
        negLL=0;
        [S,~,C]=unique(obj.stdFitIndParam);

        for i = 1:numel(S)
            rInd=find(obj.stdFitIndParam==S(i));
            s=obj.st(rInd);
            s=s(1);

            Rcmp=obj.RCmpChs(:,rInd,:);
            cmpx=obj.cmpR(:,rInd);
            nl=obj.negLLFunStd(Rcmp,cmpx,s);
            negLL=negLL + nl;
        end
    end
%- INSERT
    function F=insert_params_2(obj,F,bInit)
        % Rho
        if bInit || obj.bRhoFitAny
            F.Rho(F.indRhoU)=obj.rho;
            F.Rho=F.Rho + F.Rho' + F.eye;
        end

        if bInit || obj.bMuFitAny
            % MU
            if bInit || obj.bMuFitAny
                F.Mu(:)=obj.mustd(:);
            end

            % Crit
            if bInit || obj.bCrFitAny
                % matricize crit
                for i = 1:obj.nPass
                    F.CRB(F.CRB==i)=obj.cr(i);
                end

                F.CRB=bsxfun(@minus,F.CRB,F.Mu);
            end
        end

    end
end
end
