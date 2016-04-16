function [LOO1,LOO2] = LOO(lamda,gamma,eta)
% [LOO1,LOO2] = LOO(0.7, 0.5, 0.5)
load Mim5NN
clear MimW
load PPIM
PPIM = (PPIM>0);
load BridgeM
load NeighboringGenes

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


[idxMIM, idxG] = find(bridgeM);
cnt = 0;
t = cputime;
Nstep = [];
% leave each phenotype-gene relationship out
for i = 1 : length(idxMIM)
    d0 = zeros(Nd,1);
    d0(idxMIM(i)) = 1; % seed phenotype
    bridgeM(idxMIM(i),idxG(i)) = 0; % remove the phenotype- gene relationship

% to calculate the transition matrix from gene network to phenotype network
    idx1 = find(sum(bridgeM) > 0);    
    G2P = bridgeM; G2P(:) = 0; % to initialize the transition matrix
    for ii = 1 : length(idx1)
        G2P(:,idx1(ii)) = bridgeM(:,idx1(ii))/sum(bridgeM(:,idx1(ii)));        
    end

% to calculate the transition matrix from phenotype network to gene network
    B = bridgeM';
    idx2 = find(sum(B) > 0);
    P2G = bridgeM'; P2G(:) = 0;
    for ii = 1 : length(idx2)
        P2G(:,idx2(ii)) = B(:,idx2(ii))/sum(B(:,idx2(ii))); 
    end
    clear B idx1  idx2 ii 

%     to give the initial value to genes 
    p0=zeros(Ng,1);
    tem = P2G(:,idxMIM(i));
    train_idx = find(tem>0); % seed genes
    if ~isempty(train_idx)    
        p0(train_idx) = 1;
        p0 = p0/sum(p0);
    end
    
    [p,d,steps] = rwrH(PPIW,MimW,G2P,P2G,gamma,lamda,eta,d0,p0);
    p(train_idx)=0; % p values for seed nodes are set at zero
    
%     to find the position of the left-out gene among test genes
    test_genes = NeighboringGenes{i};   
    test_idx = find(ismember(cell2mat(genes(:,1)),test_genes));
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