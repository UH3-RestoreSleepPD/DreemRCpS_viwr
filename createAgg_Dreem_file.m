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
load(matfileL2{1},'FileMade');

procTABV = FileMade(:,["Accel_State","StreamTime","Accel_X",...
    "Accel_Y","Accel_Z"]);
clearvars("FileMade");

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
epochBLOCKs = round(height(dataEDF)/30);
channelNames = dataEDF.Properties.VariableNames;
sampFreqs = zeros(1,length(channelNames));
for si = 1:length(channelNames)
    sampFreqs(si) = length(dataEDF.(channelNames{si}){1});
end
chanFreqs = sampFreqs;
blockend = height(dataEDF);

startINDS = transpose(1:30:blockend);
endINDS = [startINDS(2:end) - 1; blockend];
allINDS = [startINDS , endINDS];



rawHpyTab = hypnoTab2;


sleepStagTab = table(sleepSTAGES,'VariableNames',{'sleepstage'});



% Save table file

cd(saveLoc)
saveName = '';
save(saveName,"sleepStagTab");


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
    keyboard
end

outTable = [procTABV2 , hypnoTab2];

% elseif summit_startIdx > dreem_startIdx
%     hypnoTab2 = hypnoTab;
%     procTABV2 = procTABV(dreem_startIdx:end,:);
% end

end