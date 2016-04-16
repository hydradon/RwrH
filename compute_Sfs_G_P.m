load PPIM;
load bridgeM;
load Mim5NN;

Np = size(MimIDs_5080, 1);  %5080
Ng = length(genes);         %8919
Sfs_G_P_2 = zeros(Np, Ng);

for i = 1 : Np
    p_p_neighbors = find(MimM(i, :) >= 0.5);
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
        %Finding Nv - Nu for phenotypes
        NvMinusNu_p = length(setdiff(g_p_neighbors, p_p_neighbors));
        %Total Nv - Nu
        NvMinusNu = NvMinusNu_g + NvMinusNu_p;
        
        Sfs_G_P_2(i, j) = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
        
        disp(['i is' num2str(i) ', j is' num2str(j)]); 
    end
end

save Sfs2.mat Sfs_G_P_2