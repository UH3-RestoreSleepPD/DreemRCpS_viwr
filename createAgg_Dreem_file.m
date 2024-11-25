function [] = createAgg_Dreem_file(dataLoc , saveLoc)

cd(dataLoc)

% Locate and load EDF file
edfFileL1 = dir('*.edf');
edfFileL2 = {edfFileL1.name};
edfFile2L = edfFileL2{1};
edfName = extractBefore(edfFile2L,'.');
curSelName = edfName;

% Locate and load proctable - extract StreamTime, Accel_State, Accel_XYZ
matfileL1 = dir('*.mat');
matfileL2 = {matfileL1.name};
load(matfileL2{1});

procTABV = CombinedProcTable(:,["Accel_State","StreamTime","Accel_X",...
    "Accel_Y","Accel_Z"]);
clearvars("CombinedProcTable");

% Locate and load hypnogram txt file
hypnoTextN = [curSelName '_hypnogram.txt'];
hypnoTab = readtable(hypnoTextN);
% hypnoTab2 = hypnoTab(1:epochBLOCKs,:);
% Clean up table
newNAMES = {'Awake','N1','N2','N3','REM','Movement'};
oldNAMES = {'SLEEP-S0','SLEEP-S1','SLEEP-S2','SLEEP-S3',...
    'SLEEP-REM','SLEEP-MT'};
sleepSTAGES = hypnoTab.SleepStage;
for oi = 1:length(oldNAMES)
    sleepSTAGES(matches(sleepSTAGES,oldNAMES{oi})) = newNAMES(oi);
end
hypnoTab.SleepStage = sleepSTAGES;
hypnoTab = renamevars(hypnoTab,["Duration_s_","Event","Time_hh_mm_ss_"],...
    ["Duration","DreemEv","hh_mm_ss"]);

[hypoProcTable , edfKeepIND] = getDREEMepochs(procTABV,hypnoTab);

[dataEDF,~] = edfread(edfFile2L);
dataEDFnTT = timetable2table(dataEDF);  % Removes the variable named 'Time'

% epochBLOCKs = round(height(dataEDF)/30);
channelNames = dataEDF.Properties.VariableNames;
sampFreqs = zeros(1,length(channelNames));
for si = 1:length(channelNames)
    sampFreqs(si) = length(dataEDF.(channelNames{si}){1});
end
chanFreqs = sampFreqs;
blockend = height(dataEDF);

startINDS = transpose(1:30:blockend-30);
endINDS = [startINDS(2:end) - 1; blockend];
allINDS = [startINDS , endINDS];

dataNTTfns = dataEDFnTT.Properties.VariableNames;
dataNTTfns2 = dataNTTfns(~matches(dataNTTfns,'Record Time'));

% Create an empty table with the specified variable names
dataEDFnTT2 = table('Size', [0, numel(dataNTTfns2)],...
    'VariableTypes', repmat({'cell'}, 1, numel(dataNTTfns2)),'VariableNames', dataNTTfns2);

for aii = 1:height(allINDS)

    startIND = allINDS(aii,1); 
    stopIND = allINDS(aii,2);
    tmpRAWtab = dataEDFnTT(startIND:stopIND,:);

    tmpTABLE = table('Size', [0, numel(dataNTTfns2)],...
    'VariableTypes', repmat({'cell'}, 1, numel(dataNTTfns2)),'VariableNames', dataNTTfns2);

    for fnI = 1:length(dataNTTfns2)
        tmpCOL = cell2mat(tmpRAWtab.(dataNTTfns2{fnI}));
        tmpTABLE.(dataNTTfns2{fnI}){1} = tmpCOL;
    end

    dataEDFnTT2 = [dataEDFnTT2 ; tmpTABLE];

end

dataEDF2 = dataEDFnTT2(logical(edfKeepIND),:);

positionNEW = zeros(height(dataEDF2),1);
for ppi = 1:height(dataEDF2)
    tmpROW = dataEDF2.Positiongram{ppi};

    if isscalar(unique(tmpROW))
        positionNEW(ppi) = unique(tmpROW);
    else
        countALL = tabulate(tmpROW);
        [~ , mostIND] = max(countALL(:,3));
        positionNEW(ppi) = countALL(mostIND,1);
    end

end

dataEDF2.Positiongram = positionNEW;

% Clean up respiration
for acI = 1:height(dataEDF2)

    respX = dataEDF2.RespirationX{acI};
    respY = dataEDF2.RespirationY{acI};
    respZ = dataEDF2.RespirationZ{acI};

    combinedTrace = sqrt(respX.^2 + respY.^2 + respZ.^2);
    dataEDF2.RespACCEL{acI} = combinedTrace;
end

dataEDF2 = removevars(dataEDF2,["RespirationX","RespirationY","RespirationZ"]);

% Clean up AccelState
for acS = 1:height(hypoProcTable)
    tmpACELs = hypoProcTable.Accel_State{acS};

    if ~isstruct(tmpACELs)
        hypoProcTable.RCPS_Accel{acS} = NaN;
        continue
    else
        tmpFnames = fieldnames(tmpACELs);
        if matches('accelData',tmpFnames)
            rawACC = tmpACELs.accelData;

            meanR = mean(rawACC);
            stdR = std(rawACC);

            zNData = (rawACC - meanR) ./ stdR;

            combTzn = sqrt(zNData(:,1).^2 +...
                zNData(:,2).^2 +...
                zNData(:,3).^2);

            combTznSm = smoothdata(combTzn,'gaussian',20);

            hypoProcTable.RCPS_Accel{acS} = combTznSm;
        else
            hypoProcTable.RCPS_Accel{acS} = NaN;
        end
    end
end

hypoProcTable = removevars(hypoProcTable,"Accel_State");

dreemEEG_HPY_ACCEL = [hypoProcTable , dataEDF2];

sampChtab = table(transpose(channelNames), transpose(sampFreqs),...
    'VariableNames',{'Channel','SamplingHz'});

% Save table file

cd(saveLoc)
saveName = [curSelName,'__HPYACCEL.mat'];
save(saveName,"dreemEEG_HPY_ACCEL","sampChtab",'-v7.3');


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB-Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outTable , edfINDEX] = getDREEMepochs(procTABV,hypnoTab)

I = ~cellfun(@isnan,  procTABV.StreamTime );

lfp_start = cellfun( @(x) datetime(x/1000,'ConvertFrom','posixtime'),...
    procTABV.StreamTime(I) ) - hours(5) - procTABV.Time(I);
lfp_start = mean(lfp_start);

dreem_start =  dateshift( lfp_start, 'start','day') + hypnoTab.hh_mm_ss(1);
offset = round(seconds(dreem_start - lfp_start) / 30);
if( offset > 0 )
    summit_startIdx = offset;
    dreem_startIdx = 1;
elseif( offset < 0 )
    summit_startIdx = 1;
    dreem_startIdx = -offset;
else
    summit_startIdx = 1;
    dreem_startIdx = 1;
end

% if summit_startIdx < dreem_startIdx
hypnoTab2 = hypnoTab(summit_startIdx+1:end,:);
procTABV2 = procTABV(dreem_startIdx:end,:);

edfINDEX = ones(height(hypnoTab),1);
edfINDEX(1:summit_startIdx) = 0;

if height(hypnoTab2) ~= height(procTABV2) 
    % lots of troubleshooting on the horizon
    if height(hypnoTab2) > height(procTABV2)

        edfINDEX = ones(height(hypnoTab),1);
        interval = seconds(30);
        lfpEPOCHs = height(procTABV2);

        % Generate the vector of times
        lfptVec = lfp_start + (0:lfpEPOCHs-1)' * interval;

        hypEPOCHS = height(hypnoTab2);
        hpytVec = dreem_start + (0:hypEPOCHS-1)' * interval;

        ismember(hpytVec,lfptVec)

        hpytVec.Format = 'yyyy-MM-dd HH:mm';
        lfptVec.Format = 'yyyy-MM-dd HH:mm';

        [~,minLOC] = min(abs(hpytVec - lfptVec(1)));

        hpytVec2 = hpytVec(minLOC:end);

        hypnoTab2 = hypnoTab2(minLOC:end,:);

        
        edfINDEX(1:minLOC) = 0;




    else
        edfINDEX = ones(height(hypnoTab),1);
        interval = seconds(30);
        lfpEPOCHs = height(procTABV2);

        % Generate the vector of times
        lfptVec = lfp_start + (0:lfpEPOCHs-1)' * interval;

        hypEPOCHS = height(hypnoTab2);
        hpytVec = dreem_start + (0:hypEPOCHS-1)' * interval;

        ismember(hpytVec,lfptVec)

        hpytVec.Format = 'yyyy-MM-dd HH:mm:ss';
        lfptVec.Format = 'yyyy-MM-dd HH:mm:ss';

        [~,minLOC] = min(abs(lfptVec - hpytVec(1)));

        minLOCa = minLOC - 1;

        procTABV2 = procTABV2(minLOC:end,:);
        hpytVec = hpytVec(1:end-1);
        hypnoTab2 = hypnoTab2(1:end-1,:);
        edfINDEX([1:minLOCa numel(edfINDEX)]) = 0;


    end
end

outTable = [procTABV2 , hypnoTab2];

% elseif summit_startIdx > dreem_startIdx
%     hypnoTab2 = hypnoTab;
%     procTABV2 = procTABV(dreem_startIdx:end,:);
% end

end