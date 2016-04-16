function cnt = abInitioSFS()
% cnt = abInitioSFS()

load CD_G_P;
load Sfs_G_P;
load PPIM;
load bridgeM;
load Mim5NN;

Ng = length(genes);
Nd = size(MimIDs_5080,1); 

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
    disp(i);
    %find the index in [1-5080] of phenotype gMim(i)
    idxD = find(ismember(MimIDs_5080,gMim(i)));
    
    idxG = find(bridgeM(idxD,:)>0); % index of phenotype ralated genes
    bridgeM(idxD,idxG) = 0; %remove all connections to genes

    %recalculate Sfs distance between phenotype idxD and genes idxG
    p_p_neighbors = find(MimM(idxD, :));
    if (MimM(idxD, idxD) == 0)
        p_p_neighbors(end + 1) = idxD; %add self-interaction if missing
    end
    p_g_neighbors = 0;
    %idxG ~ p_g_neighbors
    
    %get the Sfs distance from phenotype idxD to all genes
    FS_scores = Sfs_G_P(idxD, :);
%     CD_scores = CD_G_P(idxD, :);
    
    for j = 1 : length(idxG)
%         disp(j)
        g_p_neighbors = find(bridgeM(:, idxG(j)));
        g_g_neighbors = find(PPIM(:, idxG(j)));
        
         if (PPIM(idxG(j), idxG(j)) == 0)
             g_g_neighbors(end + 1) = idxG(j);
         end
        
        %Finding union of 1st phenotypes neigbors
        NuORNv_p = length(union(p_p_neighbors, g_p_neighbors));
        %Finding union of 1st genes neigbors
        NuORNv_g = length(union(p_g_neighbors, g_g_neighbors));
        %Total union
        NuORNv = NuORNv_p + NuORNv_g;
         
        %Replace the previous FS_score with new one %find number common phenotypes neighbors between pheno i and gene j
        NuANDNv_p = sum(ismember(p_p_neighbors, g_p_neighbors));
        
        %find number common genes neighbors between pheno i and gene j
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
        %Finding Nv - Nv for phenotypes
        NvMinusNu_p = length(setdiff(g_p_neighbors, p_p_neighbors));
        %Total Nv - Nu
        NvMinusNu = NvMinusNu_g + NvMinusNu_p;
        
        FS_score = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));;
        FS_scores(idxG(j)) = FS_score;

        %Symmetric different
        NuDeltaNv = NuMinusNv + NvMinusNu;
        
%         CD_score = NuDeltaNv/(NuANDNv + NuORNv);
%         CD_scores(idxG(j)) = CD_score;
        
    end
    
    %[p,d,steps] = rwrH(PPIW,MimW,G2P,P2G,gamma,lamda,eta,d0,p0);
    result_p = sort(FS_scores,'descend');
%     result_p = sort(CD_scores,'descend');
    if length(idxG)>1
        rankTem=[];
        for j = 1 : length(idxG)
            rankTem(j,1) = round(mean(find(result_p == FS_scores(idxG(j)))));            
%             rankTem(j,1) = round(mean(find(result_p == CD_scores(idxG(j)))));            

        end
        cnt = cnt + sum(rankTem==1);
        rank{i,1} = rankTem;  
    else
        rank{i,1} = round(mean(find(result_p == FS_scores(idxG))));
%         rank{i,1} = round(mean(find(result_p == CD_scores(idxG))));
        cnt = cnt + sum(rank{i,1}==1);
    end
    bridgeM(idxD,idxG)  = 1; %restore the gene-phenotype link
%     Nstep(i,1) = steps;
end
TTT = cputime-t
datestr(now)
