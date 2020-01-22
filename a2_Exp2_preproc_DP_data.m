% Exp3_DPanalysis
%
% last modified 13/01/2020

% clear
% close all

subs = [13,14,16:19,21:25,27:29,32:37];


cd 'D:\Behavior_3-Jitter_Control\Jitter_Control_V2.2\';
data = [pwd '\analysis\allSortedData\'];
resOut = [pwd '\results\'];

% fs = 44100;
Exp2_allDP = zeros(1,2);
% n = 20; % number of bootstraps


rhythm = [1 2];
allRTs = zeros(30,length(rhythm),length(subs));
allOutliers = zeros(79,length(subs));
for sIdx = 1:length(subs)
    
    % load results
    sub = subs(sIdx);
    %     allDataJitter(sIdx).RTs = zeros(30,2)
    load([data,'S',num2str(sub),'_run1'],'results')
    
    Exp2_allData(sIdx).results = results;
    
    
    %     figure
    %     boxplot(allDataJitter(sIdx).results(:,6))
    
    
    %% remove outliers per sub for all trials
    % remove first trial
    results = results(2:end,:);
    
    
    %define RT and corrRT
    RT = results(:,6);
    TS=results(:,3);
    
    % make no targ trials as if target after 12s
    TargPos = results(:,4);
    %     TargPos(TargPos == 0)=12;
    RT(TargPos == 0)=nan;
    
    
    corrRT = RT - TargPos;
    corrSDT = results(:,5);
    
    hits = corrSDT==3;
    temp = corrRT;
    temp(~hits)=nan;
    numHit = sum(~isnan(temp));
    corrTS = results(:,3);
    %%
    subSD = nanstd(temp);
    subMu = nanmean(temp);
    
    % determine cutoff
    subCut = subMu+(3*subSD);
    %     subCut = 3
    
    AllSc(sIdx,:) = subCut;
    jitterExp_allRT = corrRT;
    valRTnum(sIdx) = sum(~isnan(corrRT));
    
    jitterExp_allRT(jitterExp_allRT<0.1)=nan;
    toofast(sIdx) = valRTnum(sIdx)-sum(~isnan(jitterExp_allRT));
    jitterExp_allRT(jitterExp_allRT>subCut)=nan;
    tooslow(:,sIdx) = valRTnum(sIdx)-sum(~isnan(jitterExp_allRT));- toofast(sIdx)
    
    %remove outlier
    
    corrTS(jitterExp_allRT>subCut)=99;
    corrTS(jitterExp_allRT<0.1)=99;
    corrSDT(jitterExp_allRT>subCut)=0;
    corrSDT(jitterExp_allRT<0.1)=0;
    
    corrRT(jitterExp_allRT>subCut)=nan;
    corrRT(jitterExp_allRT<0.1)=nan;
    %
    
    Exp2_allData(sIdx).results(2:end,9) = corrRT;
    
    for rIdx =1:length(rhythm)
        
        
        % index rhythm condition
        trs = find(results(:,2) == rhythm(rIdx));
        % index all hits
        Allhits = find(corrSDT==3);
        AllfalseAlarm= find(corrSDT==1);
        
        % hits  in specific rhythm condition
        HitIR  = intersect(Allhits, trs);
        % FAs in specific rhythm condition
        FaIR = intersect(AllfalseAlarm,trs);
        
        
        % find numtarg in rhythm condition
        x1=  find(corrTS>0); %all TS and outlier
        x2=find(corrTS <90) ; % all 0 trials and TS
        
        allTS =intersect(x1,x2);    % intersection of x1 and x2 i.e. all TS
        numTarg = intersect(allTS,trs);
        numTarg = length(numTarg);
        
        
        %         find num nontarg in rhythm condition
        all0TS= find(corrTS==0);
        numNonTarg = intersect(all0TS,trs);
        numNonTarg = length(numNonTarg);
        %
        %         if ~(numTarg+numNonTarg+OutNum(:,sIdx)) == 40
        %             msg = 'Error. Targets, nonTargets and Outlier don''t add up Sub:';
        %             error(msg)
        %             disp(sIdx)
        %         end
        
        HitNum = numel(HitIR);
        FANum = numel(FaIR);
        
        % to avoid infinity
        % halfHit = 0.5/numTargets; old without correct for removed FAs
        halfHit = 0.5/numTarg;
        halfFA = 0.5/numNonTarg;
        
        
        
        HitProb = HitNum./numTarg;
        probFA = FANum/numNonTarg;
        
        %       Adjust only the extreme values by replacing rates of 0 with 0.5/n
        %       and rates of 1 with (n?0.5)/n where n is the number of signal or
        %       noise trials (Macmillan & Kaplan, 1985)
        
        if HitProb == 1
            HitProb = 1-halfHit;
        elseif HitProb == 0
            HitProb = halfHit;
        end
        
        if probFA == 1
            probFA = 1-halfFA;
        elseif probFA == 0
            probFA = halfFA;
        end
        
        % Calculate d-prime
        
        zHit = norminv(HitProb);
        zFA  = norminv(probFA) ;
        d(rIdx) = zHit - zFA ;
        %
        if isnan(d)
            a=1;
        end
        
        %         % bootstrapping
        %
        %         TargVec = vertcat(ones(numTarg,1),zeros(numNonTarg,1));
        %         RespVec = vertcat(ones(HitNum+FANum,1),zeros(length(TargVec)-length(HitIR)-length(FaIR),1));
        %         [t3,bootsamp]  =  bootstrp(n,@fc_dPrimeBiasBoot,TargVec,RespVec);
        %         %         btstrpData(sIdx).bootsamp(:,rIdx) = bootsamp;
        %         btstrpData(sIdx).bootstat(:,rIdx) = t3;
        %
        
        % combine all d' values in one matrix
        Exp2_allDP(sIdx,:)= d;
        
        
        % save hits/FAs
        subHits(rIdx,1)=HitNum;
        subFAs(rIdx,1)=FANum;
        
        %save targ/notarg
        subTargs(rIdx,1)= numTarg;
        subNoTargs(rIdx,1)=numNonTarg;
        
    end
    
    allHits(:,sIdx) = subHits
    allFAs(:,sIdx) = subFAs
    allTargs(:,sIdx) = subTargs
    allNoTargs(:,sIdx) = subNoTargs
    
    %     save([resOut,'S',num2str(sub),'_dPrime'],'d');
end

% save([resOut,'allDP_results'],'allDP_jittered','btstrpData');
save([resOut,'Exp2_allDP'],'Exp2_allDP');

save([resOut,'Exp2_SDTscores'],'allTargs','allNoTargs','allHits','allFAs');






