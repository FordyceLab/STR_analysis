function [squaredDiffTot] = diffSquaredCoreRepeatModelFromKD(params,koffData,KDData,weights,seqIdentifier,paramIndex)
%paramsCurr = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
%params =  [%konmax %koff,M all_koffmicro_core all_koffmicro_flank all_p_core_rel all_p_flank_rel]

[numInData,numReps]=size(koffData);

[KDModel,koffModel] = getManyCoreRepeatModel(params,seqIdentifier,paramIndex);

squaredDiffs=zeros(numInData,numReps);

for i=1:numInData
    for j=1:numReps
        squaredDiffs(i,j)=(((KDModel(i)-KDData(i,j))/KDData(i,j))^2+((koffModel(i)-koffData(i,j))/koffData(i,j))^2)*weights(i);
    end
end
    squaredDiffTot=nansum(nansum(squaredDiffs));

end



