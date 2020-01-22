%% Exp4 (aperiodic, periodic and cue experiment) 

% this script reads the raw data, removes RT (and hit) outlier per subejct and saves SDT
% and RT scores in matrices for further analyses

% last modifed 22/08/2019


% load all runs per sub, combine in matrix

% clear
% close all

subs = [1,3:4,6:12, 14,15,18:20,22:25];
NumRun = [1:6];

cd 'D:\Behavior_4-CueExp\';
outdir = [pwd '\experiment\Output\'];
resOut = [pwd '\results\'];

% preallocate
data = zeros(24,7,3,2);
Exp4_allDP = zeros(3,length(subs));
Exp4_allRT = zeros(numel(NumRun)*16,3,length(subs)); % hardcoded, 16 = number of target trials per run per condition


for sIdx = 1:length(subs)
    
    % load raw data
    sub = subs(sIdx);
    
    for iR= 1:numel(NumRun)
        load([outdir,'S',num2str(sub),'_Block',num2str(NumRun(iR)),'_results']);
        data(:,:,:,iR) = results;
    end
    
    % raw data (24 trials x 3 conditions x 6 runs) 
    RT = squeeze(data(:,5,:,:));        % reaction time relative to start of targetQuintet
    TargPos = squeeze(data(:,2,:,:));   % number of target quintet within sequence
    SDT = squeeze(data(:,6,:,:));       % hits/misses/CR/FalseAlarms (3 = hit; 1 = fa) 
    
    %% remove outliers
    
    % sort all trials by condition per column
    for idx = 1:3
        temp = squeeze(RT(:,idx,:));
        temp = temp(:);
        catRT(:,idx) = temp;
        
        temp = squeeze(SDT(:,idx,:));
        temp = temp(:);
        catSDT(:,idx) = temp;
        
        temp = squeeze(TargPos(:,idx,:));
        temp = temp(:);
        catTargPos(:,idx) = temp;
    end
    
    %to check below whether correct trials removed
%     catSDT2 = catSDT;
  
    %all trials classified as hits are used for RT analysis 
    hits = catSDT==3;
    %number of hits, without taking RT int consideration yet
    uncorrHits(:,sIdx) =  sum(hits);
    
    for idx = 1:3
        corrRT = catRT(hits(:,idx),idx);
        Exp4_allRT(1:length(corrRT),idx,sIdx) = corrRT;
        valRTnum(idx,sIdx) = length(corrRT);
    end
    
    % use these hit trials to determine RT outlier cutoff 
    catRT(~hits)=nan;
    % mean and std of reaction time across all trials per subject
    subSD = nanstd(catRT(:));
    subMu = nanmean(catRT(:));
    
    % determine cutoff
    subCut = subMu+(3*subSD);
    
    AllSc(sIdx,:) = subMu;
  
    Exp4_allRT(Exp4_allRT<0.1)=nan;
    toofast(:,sIdx) = valRTnum(:,sIdx)'-sum(~isnan(Exp4_allRT(:,:,sIdx)));
    Exp4_allRT(Exp4_allRT>subCut)=nan;
    tooslow(:,sIdx) = valRTnum(:,sIdx)'-sum(~isnan(Exp4_allRT(:,:,sIdx)));
    valTrials= squeeze(sum(~isnan(Exp4_allRT)));
    
    
    %% DP analysis cutoff 
    % use RT outlier cutoff to remove hits that were too slow or fast 
    
    catSDT(catRT>subCut)=nan;
    catSDT(catRT<0.1)=nan;
   
    
    %check whether corrct trials are sorted out 
%     sortOut(:,sIdx) = sum(isnan(catSDT2-catSDT));
%     sum(sortOut');
    
    % per condition (periodic, aperiodic, cue)
 
    h = catSDT==3;%find hit trials
    fa = catSDT ==1;%find false alarm trials
    
    
    for conIdx = 1:3
        
        %identify targets and no targets per run
        zeroes = catTargPos == 0;
        ones =  ~isnan(catTargPos) & catTargPos ~= 0;
        
        noTargs = sum(zeroes(:,conIdx)); %48
        Targs = sum(ones(:,conIdx)); %96
        
        % number of hits and FAs in condition (each element in vector is number
        % of hits/ FAs per run for a subject
        falseA = sum(fa(:,conIdx));
        hit = sum(h(:,conIdx));
        
       
        %% calculate d prime for plotting 
        
        % edge correction: 
        % adjust only the extreme values by replacing rates of 0 with 0.5/n
        % and rates of 1 with (n-0.5)/n 
        % where n is the number of signal or noise trials (Macmillan & Kaplan, 1985)
        % to avoid infinity
        
        halfHit = 0.5./Targs;
        halfFA = 0.5./noTargs;
        
        HitProb = hit./Targs;
        probFA = falseA./noTargs;

        if HitProb == 1
            HitProb = 1-halfHit;
        elseif HitProb == 0
            HitProb = halfHit;
        end

        if probFA== 1
            probFA = 1-halfFA;
        elseif probFA == 0
            probFA = halfFA;
        end
        
        zHit = norminv(HitProb);
        zFA  = norminv(probFA) ;
        d(conIdx,:) = zHit - zFA ;
        
        
        % save hits/FAs
        subHits(conIdx,1)=hit;
        subFAs(conIdx,1)=falseA;
        
        %save targ/notarg
        subTargs(conIdx,1)= Targs;
        subNoTargs(conIdx,1)=noTargs;
        
        
    end
    
    % save Hits False Alarms and number of target/ no target trials
    allHits(:,sIdx) = subHits;
    allFAs(:,sIdx) = subFAs;
    allTargs(:,sIdx) = subTargs;
    allNoTargs(:,sIdx) = subNoTargs;

    save([resOut,'SDTscores'],'allTargs','allNoTargs','allHits','allFAs');
    
    
    Exp4_allDP(:,sIdx)= d;
%     Exp4_allRT(:,:,sIdx)= catRT;
end

save([resOut,'Exp4_allDP_ConXSub'],'Exp4_allDP');
save([resOut,'Exp4_allRT_ConXSub'],'Exp4_allRT');


