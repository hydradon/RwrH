load PPIM;
% load Mim5NN;

Ng = length(genes);
% Nd = length(5080);
Sfs_G_G = zeros(Ng);
% Sfs_P_P = zeros(Nd);

 for i = 1 : (Ng - 1)
    for j = (i+1) : Ng
%         if (i < 5081) & (j < 5081)
%              Nu = find(MimM(i, :));
%             if (MimM(i, i) == 0)
%                 Nu(end + 1) = i;
%             end
% 
%             Nv = find(MimM(j, :));
%             if (MimM(j, j) == 0)
%                  Nv(end + 1) = j;
%             end
% 
%             NuANDNv = sum(ismember(Nu, Nv));
%             NuMinusNv = length(setdiff(Nu, Nv));
%             NvMinusNu = length(setdiff(Nv, Nu));
% 
%             Sfs_P_P(i, j)   = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
%             Sfs_P_P(j, i) = Sfs_P_P(i, j);
%             
%         end
        
        Nu = find(PPIM(i, :));
        if (PPIM(i, i) == 0)
            Nu(end + 1) = i;
        end
        
        Nv = find(PPIM(j, :));
        if (PPIM(j, j) == 0)
             Nv(end + 1) = j;
        end
        
        NuANDNv = sum(ismember(Nu, Nv));
        NuMinusNv = length(setdiff(Nu, Nv));
        NvMinusNu = length(setdiff(Nv, Nu));
        
        Sfs_G_G(i, j)   = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
        Sfs_G_G(j, i) = Sfs_G_G(i, j);
        
        disp(['i is' num2str(i) ', j is' num2str(j)]);
    end
end

save Sfs_G_G.mat Sfs_G_G
% save Sfs_P_P.mat S fs_P_P
