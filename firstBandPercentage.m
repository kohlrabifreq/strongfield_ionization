clear all;
latticeDepth = 8; %in recoils
pulseFreq = 0.1:0.05:6; %in kHz
firstBand = [];

for mm=1:length(pulseFreq)
    output = strongFieldAnalysisV2(latticeDepth, pulseFreq(mm));
    firstBand = [firstBand output.FirstBandPercent];
    disp(['Frequency: ', num2str(output.FreqkHz), 'kHz']);
end

plot(pulseFreq,firstBand);