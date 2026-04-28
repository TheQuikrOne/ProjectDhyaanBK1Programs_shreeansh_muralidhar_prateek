% pairLevelGamma - Per-pair Med-Con direction counts for new gamma pipeline.
addpath(genpath('C:\AIML\patagonia\nsp_project'));

savedDataFolder = fullfile(fileparts(mfilename('fullpath')), 'savedDataGamma');
pairedSubjectNameList = getPairedSubjectsBK1;
protocolsToAnalyze = {'EO1','M1'};
numPairs = size(pairedSubjectNameList,1);
numProtocols = length(protocolsToAnalyze);

kMed=NaN(numPairs,numProtocols); kCon=kMed;
velMed=kMed; velCon=kMed; volMed=kMed; volCon=kMed;
radMed=kMed; radCon=kMed; tnnMed=kMed; tnnCon=kMed; prMed=kMed; prCon=kMed;

for pairIdx=1:numPairs
    for grpIdx=1:2
        for pIdx=1:numProtocols
            fname = fullfile(savedDataFolder, ...
                [pairedSubjectNameList{pairIdx,grpIdx} '_' protocolsToAnalyze{pIdx} '_stateSpace.mat']);
            if ~exist(fname,'file'); continue; end
            r = load(fname);
            if grpIdx==1
                kMed(pairIdx,pIdx)=r.k; velMed(pairIdx,pIdx)=r.metrics.velocity;
                volMed(pairIdx,pIdx)=r.metrics.logVolume; radMed(pairIdx,pIdx)=r.metrics.radius;
                tnnMed(pairIdx,pIdx)=r.tnnD;
                if isfield(r.metrics,'PR'); prMed(pairIdx,pIdx)=r.metrics.PR; end
            else
                kCon(pairIdx,pIdx)=r.k; velCon(pairIdx,pIdx)=r.metrics.velocity;
                volCon(pairIdx,pIdx)=r.metrics.logVolume; radCon(pairIdx,pIdx)=r.metrics.radius;
                tnnCon(pairIdx,pIdx)=r.tnnD;
                if isfield(r.metrics,'PR'); prCon(pairIdx,pIdx)=r.metrics.PR; end
            end
        end
    end
end

metrics = {'k',kMed,kCon; 'Velocity',velMed,velCon; 'LogVolume',volMed,volCon; ...
           'Radius',radMed,radCon; 'PR',prMed,prCon; 'TNN_d',tnnMed,tnnCon};

for pIdx = 1:numProtocols
    fprintf('\n=== Pair-level direction counts at %s (Med - Con) ===\n', protocolsToAnalyze{pIdx});
    fprintf('%-10s %8s %8s %8s %10s %12s\n','Metric','M>C','M<C','M=C','% M>C','Binomial p');
    fprintf('%s\n', repmat('-',1,60));
    for m=1:size(metrics,1)
        M = metrics{m,2}(:,pIdx); C = metrics{m,3}(:,pIdx);
        v = ~isnan(M) & ~isnan(C); n = sum(v);
        d = M(v) - C(v);
        nGreater = sum(d > 0); nLess = sum(d < 0); nEq = sum(d == 0);
        nNonZero = nGreater + nLess;
        % Two-tailed exact binomial test against 0.5
        if nNonZero > 0
            k = min(nGreater, nLess);
            p = 2 * binocdf(k, nNonZero, 0.5);
            p = min(p, 1);
        else
            p = NaN;
        end
        fprintf('%-10s %8d %8d %8d %9.1f%% %12.4f\n', ...
            metrics{m,1}, nGreater, nLess, nEq, 100*nGreater/n, p);
    end
end
