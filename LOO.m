function [LOO1,LOO2] = LOO(lamda,gamma,eta)
% [LOO1,LOO2] = LOO(0.7, 0.5, 0.5)
load Mim5NN %Mim5NN includes MimM and MimW, each 5080 * 8919
clear MimW
load PPIM      %PPI (A_G)matrix 8919 * 8919
PPIM = (PPIM>0);
load BridgeM                    %5080 phenotypes * 8919 genes matrix
load NeighboringGenes   %1428 * 1 (gene-phenotype links)

Ng = length(genes);                 %gene: 8919 genes * 5 matrix
Nd = size(MimIDs_5080,1);   %phenotype: 5080 * 1
% to get the transition matrix for gene network
for i = 1 : Ng
    PPIW(:,i) = PPIM(:,i)/sum(PPIM(:,i));
end
%clear PPIM

% to get the transition matrix for phenotype network
for i = 1 : Nd
    MimW(:,i) = MimM(:,i)/sum(MimM(:,i));
end
%clear MimM


[idxMIM, idxG] = find(bridgeM); %1428, 1428

cnt = 0;
t = cputime;
Nstep = [];
% leave each phenotype-gene relationship out
for i = 1 : length(idxMIM)
    d0 = zeros(Nd,1);   %return a Nd * 1 all zeros matrix
    
    %----------------------------------------------------------------------
    %assigning neighboring phenotypes as seed phenotypes
    adjPhenotype = find(MimM(:, idxMIM(i))); %6  
%    for j = 1 : length(adjPhenotype)
%         d0(adjPhenotype(j)) = MimM(adjPhenotype(j), idxMIM(i));
%    end
    %----------------------------------------------------------------------
    
   
    
    d0(idxMIM(i)) = 1; % seed phenotype
    
    d0 = rwr(MimM, d0, 0.3); %gamma = 0.1, 0.3, 0.5, 0.7, 0.9
    
%     d0 = d0/sum(d0);
    bridgeM(idxMIM(i),idxG(i)) = 0; % remove the phenotype- gene relationship

% to calculate the transition matrix from gene network to phenotype network
    %sum(bridgeM) =  rows of sum of each collumn
    idx1 = find(sum(bridgeM) > 0);    %937
    
    G2P = bridgeM; G2P(:) = 0; % to initialize the transition matrix
    for ii = 1 : length(idx1)
        G2P(:,idx1(ii)) = bridgeM(:,idx1(ii))/sum(bridgeM(:,idx1(ii)));        
    end

% to calculate the transition matrix from phenotype network to gene network
    B = bridgeM';
    idx2 = find(sum(B) > 0);  %1126
    P2G = bridgeM'; P2G(:) = 0;
    for ii = 1 : length(idx2)
        P2G(:,idx2(ii)) = B(:,idx2(ii))/sum(B(:,idx2(ii))); 
    end
    clear B idx1  idx2 ii 

%     to give the initial value to genes 
    p0=zeros(Ng,1);
    tem = P2G(:,idxMIM(i))
    train_idx = find(tem>0); % seed genes
    
    if ~isempty(train_idx)    
        p0(train_idx) = 1;
        
        %----------------------------------------------------------------------
        %find neighboring genes of seed genes
%         for j = 1 : length(train_idx)
%             adjGenes = find(PPIM(:, train_idx(j)));
%     
%             for a = 1 : length(adjGenes)
%                 if (p0(adjGenes(a)) ~= 1)
%                     p0(adjGenes(a)) = PPIM(adjGenes(a), train_idx(j));
%                 end
%             end
%         end
        %----------------------------------------------------------------------
        
        
        
        %----------------------------------------------------------------------
        %assigning genes of neighboring phenotypes as seed genes
%         for j = 1 : length(adjPhenotype)
%             causeGenes = find(P2G(:, adjPhenotype(j)));
%             pSimilarity = MimM(adjPhenotype(j), idxMIM(i));
%             
%             for a = 1 : length(causeGenes)
%                 seedGeneValue = pSimilarity * P2G(causeGenes(a), adjPhenotype(j));
%                 
%                 if(p0(causeGenes(a)) < seedGeneValue)
%                     p0(causeGenes(a)) = seedGeneValue;
%                 end
%             end
%         end
        %----------------------------------------------------------------------
        
        
        p0 = p0/sum(p0);
    end
    
   
    
    [p,d,steps] = rwrH(PPIW,MimW,G2P,P2G,gamma,lamda,eta,d0,p0);
    p(train_idx)=0; % p values for seed nodes are set at zero
    
%     to find the position of the left-out gene among test genes
    test_genes = NeighboringGenes{i};   
    
    %test_idx: indices of all neighbor genes of idxG(i) gene
    test_idx = find(ismember(cell2mat(genes(:,1)),test_genes));
    
    %rank: 1428 * 2
    result_p = sort(p(test_idx),'descend');
    rank(i,1) = round(mean(find(result_p == p(idxG(i)))));
    
    bridgeM(idxMIM(i),idxG(i)) = 1; % restore briging matrix for the next loop
    
    Nstep(i,1) = steps; % to count the number steps to converge
    if rank(i,1)==1 cnt = cnt + 1; end
   if (~mod(i,10) || i > 1400)
       disp(['////////////////// ' num2str(cnt) '  in ' num2str(i) '  \\\\\\\\\\\\\\\\\\\\'])
   end  
   
  result_p2 = sort(p,'descend');
  rank(i,2) = round(mean(find(result_p2 == p(idxG(i)))));
  
%   rank(i,3) = MimIDs_5080(idxMIM(i));
%   rank(i,4) = genes{idxG(i),1}; % hprd id
end
TTT = cputime-t;
LOO1 = sum(rank(:,1)==1);
LOO2 = sum(rank(:,2)==1);
% save Results_PPI rank Nstep TTT cutoff lamda gamma Results