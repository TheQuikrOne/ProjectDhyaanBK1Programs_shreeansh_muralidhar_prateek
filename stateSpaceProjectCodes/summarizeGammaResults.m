% summarizeGammaResults - Quick aggregation of new gamma results for comparison.
addpath(genpath('C:\AIML\patagonia\nsp_project'));

savedDataFolder = fullfile(fileparts(mfilename('fullpath')), 'savedDataGamma');
pairedSubjectNameList = getPairedSubjectsBK1;
protocolsToAnalyze = {'EO1','M1'};
numPairs = size(pairedSubjectNameList,1);
numProtocols = length(protocolsToAnalyze);

K_std = NaN;
kMed=NaN(numPairs,numProtocols); kCon=kMed;
velMed=kMed; velCon=kMed; volMed=kMed; volCon=kMed;
radMed=kMed; radCon=kMed; tnnMed=kMed; tnnCon=kMed;
prMed=kMed; prCon=kMed; varMed=kMed; varCon=kMed;

for pairIdx=1:numPairs
    for grpIdx=1:2
        for pIdx=1:numProtocols
            fname = fullfile(savedDataFolder, ...
                [pairedSubjectNameList{pairIdx,grpIdx} '_' protocolsToAnalyze{pIdx} '_stateSpace.mat']);
            if ~exist(fname,'file'); continue; end
            r = load(fname);
            if isnan(K_std) && isfield(r,'K_std'); K_std = r.K_std; end
            if grpIdx==1
                kMed(pairIdx,pIdx)=r.k; velMed(pairIdx,pIdx)=r.metrics.velocity;
                volMed(pairIdx,pIdx)=r.metrics.logVolume; radMed(pairIdx,pIdx)=r.metrics.radius;
                tnnMed(pairIdx,pIdx)=r.tnnD;
                if isfield(r.metrics,'PR'); prMed(pairIdx,pIdx)=r.metrics.PR; end
                if isfield(r.metrics,'totalVariance'); varMed(pairIdx,pIdx)=r.metrics.totalVariance; end
            else
                kCon(pairIdx,pIdx)=r.k; velCon(pairIdx,pIdx)=r.metrics.velocity;
                volCon(pairIdx,pIdx)=r.metrics.logVolume; radCon(pairIdx,pIdx)=r.metrics.radius;
                tnnCon(pairIdx,pIdx)=r.tnnD;
                if isfield(r.metrics,'PR'); prCon(pairIdx,pIdx)=r.metrics.PR; end
                if isfield(r.metrics,'totalVariance'); varCon(pairIdx,pIdx)=r.metrics.totalVariance; end
            end
        end
    end
end

fprintf('\n=== GAMMA (new pipeline) K_std = %d ===\n\n', K_std);
fprintf('%-12s %-6s %16s %16s %10s\n','Metric','Proto','Med (mean+/-SEM)','Con (mean+/-SEM)','p-val');
fprintf('%s\n', repmat('-',1,70));

metrics = {'k',kMed,kCon; 'Velocity',velMed,velCon; 'LogVol',volMed,volCon; ...
           'Radius',radMed,radCon; 'PR',prMed,prCon; 'TNN_d',tnnMed,tnnCon};

for m=1:size(metrics,1)
    M = metrics{m,2}; C = metrics{m,3};
    for pIdx=1:numProtocols
        v = ~isnan(M(:,pIdx)) & ~isnan(C(:,pIdx)); n=sum(v);
        if n<5; continue; end
        p = signrank(M(v,pIdx), C(v,pIdx));
        fprintf('%-12s %-6s %7.2f +/- %5.2f %7.2f +/- %5.2f %10.4f\n', ...
            metrics{m,1}, protocolsToAnalyze{pIdx}, ...
            nanmean(M(:,pIdx)), nanstd(M(:,pIdx))/sqrt(n), ...
            nanmean(C(:,pIdx)), nanstd(C(:,pIdx))/sqrt(n), p);
    end
end

fprintf('\n=== Within-group EO1 vs M1 ===\n');
fprintf('%-12s %-6s %16s %16s %10s\n','Metric','Group','EO1 (mean+/-SEM)','M1 (mean+/-SEM)','p-val');
fprintf('%s\n', repmat('-',1,70));
for m=1:size(metrics,1)
    for grpIdx=1:2
        if grpIdx==1; D = metrics{m,2}; gl='Med'; else; D = metrics{m,3}; gl='Con'; end
        v = ~isnan(D(:,1)) & ~isnan(D(:,2)); n=sum(v);
        if n<5; continue; end
        p = signrank(D(v,1), D(v,2));
        fprintf('%-12s %-6s %7.2f +/- %5.2f %7.2f +/- %5.2f %10.4f\n', ...
            metrics{m,1}, gl, ...
            nanmean(D(:,1)), nanstd(D(:,1))/sqrt(n), ...
            nanmean(D(:,2)), nanstd(D(:,2))/sqrt(n), p);
    end
end

fprintf('\n=== ANCOVA: LogVolume ~ Group + TotalVariance ===\n');
for pIdx=1:numProtocols
    v = ~isnan(volMed(:,pIdx)) & ~isnan(volCon(:,pIdx)) & ...
        ~isnan(varMed(:,pIdx)) & ~isnan(varCon(:,pIdx));
    n=sum(v);
    if n<6; fprintf('%s: n=%d, skipped\n', protocolsToAnalyze{pIdx}, n); continue; end
    volV = [volMed(v,pIdx); volCon(v,pIdx)];
    grpV = [ones(n,1); zeros(n,1)];
    varV = [varMed(v,pIdx); varCon(v,pIdx)];
    tbl = table(volV,grpV,varV,'VariableNames',{'LogVolume','Group','TotalVariance'});
    lm = fitlm(tbl,'LogVolume ~ Group + TotalVariance');
    fprintf('\n  Protocol %s (n=%d per group):\n', protocolsToAnalyze{pIdx}, n);
    disp(lm.Coefficients);
    fprintf('  R^2 = %.3f\n', lm.Rsquared.Ordinary);
end
