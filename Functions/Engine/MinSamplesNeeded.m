function numberSamples=MinSamplesNeeded(scoreData,scoreThreshold,minConf)

numberReps=size(scoreData,2);
numberRepsLessThanThreshold=sum(scoreData<scoreThreshold,2);

numberSamples=find((numberRepsLessThanThreshold/numberReps)>minConf,1);
if(isempty(numberSamples))
    numberSamples=NaN;
end

end