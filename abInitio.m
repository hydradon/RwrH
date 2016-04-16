function cnt = abInitio(lamda,gamma,eta)
% cnt = abInitio(0.5,0.7,0.5)
load Mim5NN
% clear MimM
load PPIM
PPIM = (PPIM>0);
load BridgeM

Ng = length(genes);
Nd = size(MimIDs_5080,1); 

% to get the transition matrix for gene network
for i = 1 : Ng
    PPIW(:,i) = PPIM(:,i)/sum(PPIM(:,i));
end
clear PPIM

% to get the transition matrix for phenotype network
for i = 1 : Nd
    MimW(:,i) = MimM(:,i)/sum(MimM(:,i));
end
clear MimM

[idxMIM, idxPPI] = find(bridgeM);
gMim = unique(MimIDs_5080(idxMIM)); %MimIDs having gene annotation

cnt = 0;
t = cputime;
Nstep = [];
p0=zeros(Ng,1);
for i = 1 : length(gMim)
    idxD = find(ismember(MimIDs_5080,gMim(i)));
    idxG = find(bridgeM(idxD,:)>0); % index of phenotype ralated genes
    bridgeM(idxD,idxG) = 0;
    [G2P,P2G] = getBridgeM(bridgeM);
    d0 = zeros(Nd,1); d0(idxD) = 1; % seed phenotype
    [p,d,steps] = rwrH(PPIW,MimW,G2P,P2G,gamma,lamda,eta,d0,p0);
    result_p = sort(p,'descend');
    if length(idxG)>1
        rankTem=[];
        for j = 1 : length(idxG)
            rankTem(j,1) = round(mean(find(result_p == p(idxG(j)))));            
        end
        cnt = cnt + sum(rankTem==1);
        rank{i,1} = rankTem;  
    else
        rank{i,1} = round(mean(find(result_p == p(idxG))));
        cnt = cnt + sum(rank{i,1}==1);
    end
    bridgeM(idxD,idxG)  = 1;
    Nstep(i,1) = steps;
end
TTT = cputime-t;
datestr(now)
