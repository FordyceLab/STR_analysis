function [KDModel,koffModel] = getManyCoreRepeatModelKDmicro(params,weights,seqIdentifier,paramIndex)

%params =  [%konmax %koff,M all_koffmicro_core all_koffmicro_flank all_p_core_rel all_p_flank_rel]

[numInData]=size(seqIdentifier,1);

KDModel=NaN*zeros(numInData,1);
koffModel=NaN*zeros(numInData,1);
for i=1:numInData
    if weights(i)>0
        paramsCurr=zeros(1,6);
        paramsCurr(1)=params(1); %konmax
        paramsCurr(2)=params(2);%koffM
        coreIndex = seqIdentifier(i,1);
        if coreIndex==0
            paramsCurr(3)=1;%Shouldn matter what you set here
            paramsCurr(5)=0;
        else
            paramsCurr(3)=params(paramIndex(1)-1+coreIndex);% koffmicro_core
            paramsCurr(5)=params(paramIndex(3)-1+coreIndex);%p_core_rel
        end
        
        flankIndex = seqIdentifier(i,2);
        if flankIndex==0
            paramsCurr(4)=1;%Shouldn matter what you set here
            paramsCurr(6)=0;
        else
            paramsCurr(4)=params(paramIndex(2)-1+flankIndex);% koffmicro_core
            paramsCurr(6)=params(paramIndex(4)-1+flankIndex);%p_core_rel
        end
        if paramsCurr(6)>0||paramsCurr(5)>0
            [KDModel(i),koffModel(i)] = getRatesCoreRepeatFastConvertParamsKDmicro(paramsCurr,1);
        else
            KDModel(i)=NaN;
            koffModel(i)=NaN;
        end
    end
end
end

