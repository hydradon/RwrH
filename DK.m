function [LOO1,LOO2] = DK(beta)
%[LOO1, LOO2] = DK(0.2)
load Mim5NN
load PPIM
load BridgeM
load NeighboringGenes
% load MimAdj

Ng = length(genes);
Nd = size(MimIDs_5080, 1);

%get adjacency matrix of MimM
MimM_ADJ = zeros(Nd);
for i = 1 : Nd
    MimM_ADJ(i, :) = MimM(i, :) >= 0.2; %only pick neigbours with similarity score >= 0.2
end

[idxMIM, idxG] = find(bridgeM); %1428, 1428
% leave each phenotype-gene relationship out
for i = 1 : length(idxMIM)
    bridgeM(idxMIM(i),idxG(i)) = 0; % remove the phenotype- gene relationship
    
    huge = [PPIM, bridgeM'; bridgeM, MimM_ADJ]; %A
    D = zeros(Ng + Nd);
    
    %calculate adjacency matrix of total matrix
    for j = 1 : (Ng+Nd)
        D(j, j) = size(find(huge(:, j)), 1);
    end

    L = D - huge;
%     beta = 0.2;
    K = expm(-beta * L);
    
    %other casautive genes of idxMIM(i)
    other_genes = find(bridgeM(idxMIM(i), :));
    
    %get the probability vector of node idxMIM(i) in K
    prob_vector = K(:, idxMIM(i) + Ng);
    %get the first 8919 probabilities (8919 genes)
    scores = prob_vector(1:Ng);
    
    for x = 1 : length(other_genes)
        genes_prob_vector = K(:, other_genes(x));
        genes_scores = genes_prob_vector(1:Ng);
        scores = scores + genes_prob_vector;
    end
    
    test_genes = NeighboringGenes{i};
    
    %test_idx: indices of all neighbor genes of idxG(i) gene
    test_idx = find(ismember(cell2mat(genes(:,1)),test_genes));
    result_p = sort(p(test_idx),'descend');
    rank(i,1) = round(mean(find(result_p == p(idxG(i)))));
    
    if rank(i,1)==1 cnt = cnt + 1; end
    if (~mod(i,10) || i > 1400)
        disp(['////////////////// ' num2str(cnt) '  in ' num2str(i) '  \\\\\\\\\\\\\\\\\\\\'])
    end
    
    result_p2 = sort(scores, 'descend');
    rank(i,2) = round(mean(find(result_p2 == p(idxG(i)))));
    
    bridgeM(idxMIM(i),idxG(i)) = 1; % restore briging matrix for the next loop
end

LOO1 = sum(rank(:,1)==1);
LOO2 = sum(rank(:,2)==1);
