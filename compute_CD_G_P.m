load PPIM;
load bridgeM;
load Mim5NN;

Np = size(MimIDs_5080, 1);
Ng = length(genes);
CD_G_P = zeros(Np, Ng);

for i = 1 : Np
    p_p_neighbors = find(MimM(i, :));
    if (MimM(i, i) == 0)
        p_p_neighbors(end + 1) = i;
    end
    p_g_neighbors = find(bridgeM(i, :));
    
    for j = 1 : Ng
        g_p_neighbors = find(bridgeM(:, j));
        g_g_neighbors = find(PPIM(:, j));
        if (PPIM(j, j) == 0)
             g_g_neighbors(end + 1) = j;
        end
        
        %Finding union of 1st phenotypes neigbors
        NuORNv_p = length(union(p_p_neighbors, g_p_neighbors));
        %Finding union of 1st genes neigbors
        NuORNv_g = length(union(p_g_neighbors, g_g_neighbors));
        %Total union
        NuORNv = NuORNv_p + NuORNv_g;
        
        %find number common phenotypes neighbors between pheno i and gene j
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
        
        %Symmetric different
        NuDeltaNv = NuMinusNv + NvMinusNu;
        
%         Sfs_G_P(i, j) = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
        CD_G_P(i, j) = NuDeltaNv/(NuANDNv + NuORNv);
        
        disp(['i is' num2str(i) ', j is' num2str(j)]); 
    end
end

save CD_G_P.mat CD_G_P