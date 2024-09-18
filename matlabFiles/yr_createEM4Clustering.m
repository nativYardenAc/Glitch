function [eyeXVec,eyeYVec,MSduringVecs]=yr_createEM4Clustering()

%%% variables to change
session_name='gandalf_270618b';
cond_num=1;
mainRoot='D:\Yarden\yarden project folder\origin';
firstFrame=30;
lastFrame=140;

% mainRoot='D:\Yarden\yarden practicum files\gal';
cortexFileRoot=[mainRoot filesep 'raw_data\gandalf left\2018June27\gan_2018June27_b.1'];
calibrationFileRoot=[mainRoot filesep 'raw_data\gandalf left\2018June27\gan_2018June27_caleye.1'];

vsdfileRoot=[mainRoot filesep 'preprocessed_VSDdata' filesep session_name];
synchronyFilePath=[mainRoot filesep 'cortex-cam synched lists' filesep session_name '.xlsx'];

%%%end of definitions
firstFrInMs=(firstFrame-27).*10;
lastFrInMs=(lastFrame-27).*10;

if contains(session_name,'boromir')
    syncMat=yr_autoSyncML2MatFiles(synchronyFilePath,vsdfileRoot,cond_num);
    firstIdxSync=2;
else
    syncMat=yr_autoSyncCortex2MatFiles(synchronyFilePath,vsdfileRoot,cond_num);
    firstIdxSync=1;
end
cortex_trials=cell2mat(syncMat(firstIdxSync:end,1));
trial_nums=cell2mat(syncMat(firstIdxSync:end,4));
vsdfileRootNoisy=[mainRoot filesep 'preprocessed_VSDdata' filesep session_name filesep 'noisyfiles'];
if contains(session_name,'boromir')
    syncMat=yr_autoSyncML2MatFiles(synchronyFilePath,vsdfileRootNoisy,cond_num);
    firstIdxSync=2;
else
    syncMat=yr_autoSyncCortex2MatFiles(synchronyFilePath,vsdfileRootNoisy,cond_num);
    firstIdxSync=1;
end
cortex_trials=[cortex_trials; cell2mat(syncMat(firstIdxSync:end,1))];
trial_nums=[trial_nums; cell2mat(syncMat(firstIdxSync:end,4))];
[trial_nums_sorted,sorted_idx]=sort(trial_nums);
cortex_trials_sorted=cortex_trials(sorted_idx);
numOfTrials=trial_nums_sorted(end);
msAmpThreshold=1.5;

%eye movements params
if contains(session_name,'legolas')
    monkeySessionMetaFile.monkeyName='legolas';
    monkeySessionMetaFile.sessionName=session_name;
    monkeySessionMetaFile.engbretThreshold=7.6;
    monkeySessionMetaFile.engbertMinDur=12;
    monkeySessionMetaFile.rejectGlitch=0;
    monkeySessionMetaFile.rejectFollowers=0;
    monkeySessionMetaFile.smoothEM=0;
    monkeySessionMetaFile.smoothEM=25;
    monkeySessionMetaFile.fineTuning='hafed';
    monkeySessionMetaFile.velThreshold=8;
    monkeySessionMetaFile.ampMethod='final';
    monkeySessionMetaFile.angleMethod='vecAverage';
    monkeySessionMetaFile.msAmpThreshold=1;
    monkeySessionMetaFile.maxFrame=130;
    [MSduringVSD,mainTimeMat]=yr_msDuringVSD(cortexFileRoot,calibrationFileRoot,monkeySessionMetaFile,0);
else
    if contains(session_name,'boromir')
        monkeySessionMetaFile.monkeyName='boromir';
        monkeySessionMetaFile.sessionName=session_name;
        monkeySessionMetaFile.engbretThreshold=3.5;
        monkeySessionMetaFile.engbertMinDur=7;
        monkeySessionMetaFile.rejectGlitch=1;
        monkeySessionMetaFile.rejectInconsistent=1;
        monkeySessionMetaFile.followersMethod='reject';
        monkeySessionMetaFile.smoothEM=25;
        monkeySessionMetaFile.subSample=2;
        monkeySessionMetaFile.fineTuning='accBaseline';
        monkeySessionMetaFile.velThreshold=3;
        monkeySessionMetaFile.accThresholdBegin=2;
        monkeySessionMetaFile.accThresholdEnd=2;
        monkeySessionMetaFile.ampMethod='max';
        monkeySessionMetaFile.angleMethod='vecAverage';
        monkeySessionMetaFile.msAmpThreshold=1.5;
        monkeySessionMetaFile.maxFrame=100;
        [MSduringVSD,mainTimeMat]=yr_msDuringVSD_ml(mlFileRoot,monkeySessionMetaFile,0);
    else
        monkeySessionMetaFile.monkeyName='gandalf';
        monkeySessionMetaFile.sessionName=session_name;
        monkeySessionMetaFile.engbretThreshold=4;
        monkeySessionMetaFile.engbertMinDur=12;
        monkeySessionMetaFile.rejectGlitch=0;
        monkeySessionMetaFile.rejectInconsistent=0;
        monkeySessionMetaFile.followersMethod='ignore';
        monkeySessionMetaFile.smoothEM=25;
        monkeySessionMetaFile.fineTuning='accBaseline';
        monkeySessionMetaFile.velThreshold=3;
        monkeySessionMetaFile.accThresholdBegin=2;
        monkeySessionMetaFile.accThresholdEnd=2;
        monkeySessionMetaFile.ampMethod='max';
        monkeySessionMetaFile.angleMethod='final';
        monkeySessionMetaFile.msAmpThreshold=1.5;
        monkeySessionMetaFile.maxFrame=100;
        
        sampleRate=4;
        
        [MSduringVSD,mainTimeMat]=yr_msDuringVSD(cortexFileRoot,calibrationFileRoot,monkeySessionMetaFile,0);
    end
end

[eyeX,eyeY,time_arr,event_arr,header]=yr_calibrateCortexData(cortexFileRoot,calibrationFileRoot); 
eyeXVec=zeros((lastFrInMs-firstFrInMs)./sampleRate,numOfTrials);
eyeYVec=zeros((lastFrInMs-firstFrInMs)./sampleRate,numOfTrials);
for trial_id=1:numOfTrials
    cortex_id=cortex_trials_sorted(trial_id);
    trialFr27=cell2mat(mainTimeMat(5,cortex_id+1));
    emStartRecording=cell2mat(mainTimeMat(2,cortex_id+1));
    startEManalysis=floor((trialFr27-emStartRecording+firstFrInMs)./sampleRate);
    endEManalysis=floor((trialFr27-emStartRecording+lastFrInMs)./sampleRate);
    vec2addX=eyeX(startEManalysis:endEManalysis-1,cortex_id);
    vec2addY=eyeY(startEManalysis:endEManalysis-1,cortex_id);
    width=getfield(monkeySessionMetaFile,'smoothEM');
    if width>0
        window=round(width./sampleRate);
        if rem(window,2)==0
            window=window+1;
        end
        vecX=sgolayfilt(vec2addX,3,window);
        vecY=sgolayfilt(vec2addY,3,window);
    end
    eyeXVec(:,trial_id)=vec2addX;
    eyeYVec(:,trial_id)=vec2addY;
end

MSduringVSD_ofCond=MSduringVSD(cortex_trials_sorted,1); %msTimes in the cortex clock (ms)
MSduringVecs={"trialNum","onsetIdx","offsetIdx","amplitude (DVA)","direction (deg)","max velocity"};
ms2analyze_id=1;
%keep only MSs between frames
for trial_id=1:numOfTrials
    cortex_id=cortex_trials_sorted(trial_id);
    trialFr27=cell2mat(mainTimeMat(5,cortex_id+1));
    timeOnset=cell2mat(mainTimeMat(3,cortex_id+1))./sampleRate;
    MSsOfTrial=cell2mat(MSduringVSD_ofCond(trial_id));
    if ~isempty(MSsOfTrial)
        for ms_id=1:size(MSsOfTrial,1)
            msOnset=MSsOfTrial(ms_id,1);
            if (msOnset>=trialFr27+firstFrInMs)&&(msOnset<=trialFr27+lastFrInMs)
                MSduringVecs(ms2analyze_id+1,1)={trial_id};
                MSduringVecs(ms2analyze_id+1,2)={round((msOnset-trialFr27-firstFrInMs)./sampleRate)-1};
                MSduringVecs(ms2analyze_id+1,3)={round((msOnset-trialFr27-firstFrInMs)./sampleRate)-1+(MSsOfTrial(ms_id,2)-MSsOfTrial(ms_id,1))./sampleRate};
                MSduringVecs(ms2analyze_id+1,4)={MSsOfTrial(ms_id,3)};
                MSduringVecs(ms2analyze_id+1,5)={MSsOfTrial(ms_id,4)};
                MSduringVecs(ms2analyze_id+1,6)={MSsOfTrial(ms_id,5)};
                ms2analyze_id=ms2analyze_id+1;
            end
        end
    end
end

% for ms_id=1:size(MSduringVecs,1)-1
%     trial_id=cell2mat(MSduringVecs(ms_id+1,1));
%     cortex_id=cortex_trials_sorted(trial_id);
%     figure;
%     plot(eyeXVec(:,trial_id));
%     hold on;
%     plot(eyeYVec(:,trial_id));
%     hold on;
%     firstIdx=max(1,cell2mat(MSduringVecs(1+ms_id,2)));
%     lastIdx=min(cell2mat(MSduringVecs(1+ms_id,3)),size(eyeXVec,1));
%     plot([firstIdx:lastIdx],eyeXVec(firstIdx:lastIdx,trial_id),'r','LineWidth',2);
%     hold on;
%     plot([firstIdx:lastIdx],eyeYVec(firstIdx:lastIdx,trial_id),'r','LineWidth',2);
%     title(['cortex trial num: ' num2str(cortex_id)]);
% end

a=1;