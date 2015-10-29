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

%return unique array in sorted order
%gMim: all phenotypes that have at least one interaction to a gene
gMim = unique(MimIDs_5080(idxMIM)); %MimIDs having gene annotation
%length(gMim) = 1126 

cnt = 0;
t = cputime;
Nstep = [];
p0=zeros(Ng,1);
for i = 1 : length(gMim)
    
    %find the index in [1-5080] of phenotype gMim(i)
    idxD = find(ismember(MimIDs_5080,gMim(i)));
    
    idxG = find(bridgeM(idxD,:)>0); % index of phenotype ralated genes
    bridgeM(idxD,idxG) = 0; %remove all connections to genes
    [G2P,P2G] = getBridgeM(bridgeM);
    
    d0 = zeros(Nd,1); d0(idxD) = 1; % seed phenotype
    
%     d0 = rwr(MimM, d0, 0.3);
    
    %----------------------------------------------------------------------
    %assigning seed phenotypes
%     adjPhenotype = find(MimW(:,idxD)); 
    
%     idxMax = find(max(MimW(adjPhenotype, idxD)));
    
%     for j = 1 : length(adjPhenotype)
%         d0(adjPhenotype(j)) = MimW(adjPhenotype(j), idxD);
%             d0(adjPhenotype(j)) = 1;
%     end
%     d0(idxMax) = 1;
    %----------------------------------------------------------------------
    
    d0 = d0/sum(d0); %normalising Phenotype vector
    
%     p0=zeros(Ng,1);
    %----------------------------------------------------------------------
    %assigning causal genes of neighboring phenotypes as seeds
%      for j = 1 : length(adjPhenotype)
%             causeGenes = find(P2G(:, adjPhenotype(j)));
%             pSimilarity = MimW(adjPhenotype(j), idxD);
%             
%             for a = 1 : length(causeGenes)
%                 seedGeneValue = pSimilarity * P2G(causeGenes(a), adjPhenotype(j));
%                 
% %                 if(p0(causeGenes(a)) < seedGeneValue)
% %                     p0(causeGenes(a)) = seedGeneValue;
% %                 end
%                 p0(causeGenes(a)) = 1;
%             end
%      end
     
%      p0 = p0/sum(p0);
    %----------------------------------------------------------------------
    
    
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
TTT = cputime-t
datestr(now)
