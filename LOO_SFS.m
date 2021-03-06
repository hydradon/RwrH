function [LOO1,LOO2] = LOO_SFS()
% [LOO1,LOO2] = LOO_SFS()
load Mim5NN %Mim5NN includes MimM and MimW, each 5080 * 8919
load PPIM      %PPI (A_G)matrix 8919 * 8919
PPIM = (PPIM>0);
load bridgeM                    %5080 phenotypes * 8919 genes matrix
load NeighboringGenes   %1428 * 1 (gene-phenotype links)
load Sfs_G_P;
% load CD_G_P;

Ng = length(genes);                 %gene: 8919 genes * 5 matrix
Nd = size(MimIDs_5080,1);   %phenotype: 5080 * 1

[idxMIM, idxG] = find(bridgeM); %1428, 1428

cnt = 0;
t = cputime;
Nstep = [];
% leave each phenotype-gene relationship out
for i = 1 : length(idxMIM)
    bridgeM(idxMIM(i),idxG(i)) = 0; % remove the phenotype- gene relationship

    %recalculate Sfs distance between phenotype idxMIM(i) and genes idxG(i)
    p_p_neighbors = find(MimM(idxMIM(i), :) >= 0.0); %neighbors filtering: only choose those with >0.5 similarity scores
    if (MimM(idxMIM(i), idxMIM(i)) == 0)
        p_p_neighbors(end + 1) = idxMIM(i); %add self-interaction if missing
    end
    
    p_g_neighbors = find(bridgeM(idxMIM(i), : ));
    
    %get the Sfs distance from phenotype idxMIM(i) to all genes
    FS_scores = Sfs_G_P(idxMIM(i), :);
    %CD_scores = CD_G_P(idxMIM(i), :);
    g_p_neighbors = find(bridgeM(:, idxG(i)));
    g_g_neighbors = find(PPIM(:, idxG(i)));
    if (PPIM(idxG(i), idxG(i)) == 0)
        g_g_neighbors(end + 1) = idxG(i); %add self-interaction if missing
    end
    
    %Finding union of 1st level phenotypes neigbors
    NuORNv_p = length(union(p_p_neighbors, g_p_neighbors));
    %Finding union of 1st genes neigbors
    NuORNv_g = length(union(p_g_neighbors, g_g_neighbors));
    %Total union
    NuORNv = NuORNv_p + NuORNv_g;
    
    %find number common phenotypes neighbors between pheno idxMIM(i) and
    %gene idxG(i)
    NuANDNv_p = sum(ismember(p_p_neighbors, g_p_neighbors));    
    %find number common genes neighbors between pheno idxMIM(i) and
    %gene idxG(i)
    NuANDNv_g = sum(ismember(g_g_neighbors, p_g_neighbors));  
    %total number of common neighbors
    NuANDNv = NuANDNv_p + NuANDNv_g;
        
    %Finding Nu - Nv for genes
    NuMinusNv_g = length(setdiff(p_g_neighbors, g_g_neighbors));
    %Finding Nu - Nv for phenotypes
    NuMinusNv_p = length(setdiff(p_p_neighbors, g_p_neighbors));
    %Total Nu - Nv
    NuMinusNv = NuMinusNv_g + NuMinusNv_p;
        
    %Finding Nv - Nu for genes
    NvMinusNu_g = length(setdiff(g_g_neighbors, p_g_neighbors));
    %Finding Nv - Nu for phenotypes
    NvMinusNu_p = length(setdiff(g_p_neighbors, p_p_neighbors));
    %Total Nv - Nu
    NvMinusNu = NvMinusNu_g + NvMinusNu_p;
     
    %Symmetric different
    NuDeltaNv = NuMinusNv + NvMinusNu;
    FS_score = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));;
%     CD_score = NuDeltaNv/(NuANDNv + NuORNv);
    FS_scores(idxG(i)) = FS_score;
%     CD_scores(idxG(i)) = CD_score;

    %recalculate Sfs distance for all level-2 neighbors that use
    %idxMIM(i)-idxG(i) link
    for x = 1 : length(g_g_neighbors)
        l2_g_p_neighbors = find(bridgeM(:, g_g_neighbors(x)));
        l2_g_g_neighbors = find(PPIM(:, g_g_neighbors(x)));
        if (PPIM(g_g_neighbors(x), g_g_neighbors(x)) == 0)
            l2_g_g_neighbors(end + 1) = g_g_neighbors(x); %add self-interaction if missing
        end
        
        %use same p_p_neighbors and p_g_neighbors
        %find number common phenotypes neighbors between pheno idxMIM(i) and
        %gene g_g_neighbors(x)
        l2_NuANDNv_p = sum(ismember(p_p_neighbors, l2_g_p_neighbors));    
        %find number common genes neighbors between pheno idxMIM(i) and
        %gene g_g_neighbors(x)
        l2_NuANDNv_g = sum(ismember(l2_g_g_neighbors, p_g_neighbors));  
        %total number of common neighbors
        l2_NuANDNv = l2_NuANDNv_p + l2_NuANDNv_g;
        
        %Finding Nu - Nv for genes
        l2_NuMinusNv_g = length(setdiff(p_g_neighbors, l2_g_g_neighbors));
        %Finding Nu - Nv for phenotypes
        l2_NuMinusNv_p = length(setdiff(p_p_neighbors, l2_g_p_neighbors));
        %Total Nu - Nv
        l2_NuMinusNv = l2_NuMinusNv_g + l2_NuMinusNv_p;

        %Finding Nv - Nu for genes
        l2_NvMinusNu_g = length(setdiff(l2_g_g_neighbors, p_g_neighbors));
        %Finding Nv - Nu for phenotypes
        l2_NvMinusNu_p = length(setdiff(l2_g_p_neighbors, p_p_neighbors));
        %Total Nv - Nu
        l2_NvMinusNu = l2_NvMinusNu_g + l2_NvMinusNu_p;
        
        l2_FS_score = ((2*l2_NuANDNv)/(l2_NuMinusNv + 2*l2_NuANDNv)) * ((2*l2_NuANDNv)/(l2_NvMinusNu + 2*l2_NuANDNv));;
        FS_scores(g_g_neighbors(x)) = l2_FS_score;
    end
    %----------------------------------------------------------------

    for y = 1 : length(p_g_neighbors)
        FS_scores(p_g_neighbors(y)) = 0; %prevent remaining phenotype-related genes from being ranked top
%     CD_scores(p_g_neighbors(y)) = 0; %being ranked top
    end
    
    %[p,d,steps] = rwrH(PPIW,MimW,G2P,P2G,gamma,lamda,eta,d0,p0);
    
%     to find the position of the left-out gene among test genes
    test_genes = NeighboringGenes{i};   
    
    %test_idx: indices of all neighbor genes of idxG(i) gene
    test_idx = find(ismember(cell2mat(genes(:,1)),test_genes));
    
    %rank: 1428 * 2
    result_p = sort(FS_scores(test_idx),'descend');
    rank(i,1) = round(mean(find(result_p == FS_scores(idxG(i)))));
    
%     result_p = sort(CD_scores(test_idx),'descend');
%     rank(i,1) = round(mean(find(result_p == CD_scores(idxG(i)))));
    
    bridgeM(idxMIM(i),idxG(i)) = 1; % restore briging matrix for the next loop
    
%     Nstep(i,1) = steps; % to count the number steps to converge
    if rank(i,1)==1 cnt = cnt + 1; end
   if (~mod(i,10) || i > 1400)
       disp(['////////////////// ' num2str(cnt) '  in ' num2str(i) '  \\\\\\\\\\\\\\\\\\\\'])
   end  
   
  result_p2 = sort(FS_scores,'descend');
  rank(i,2) = round(mean(find(result_p2 == FS_scores(idxG(i)))));
  
%    result_p2 = sort(CD_scores,'descend');
%   rank(i,2) = round(mean(find(result_p2 == CD_scores(idxG(i)))));
  
%   rank(i,3) = MimIDs_5080(idxMIM(i));
%   rank(i,4) = genes{idxG(i),1}; % hprd id
end
TTT = cputime-t
LOO1 = sum(rank(:,1)==1);
LOO2 = sum(rank(:,2)==1);
% save Results_PPI rank Nstep TTT cutoff lamda gamma Results