load PPIM;
Ng = length(genes);
Sfs = zeros(Ng);

 for i = 1 : (Ng - 1)
    for j = (i+1) : Ng
        Nu = find(PPIM(i, :));
        Nv = find(PPIM(j, :));
        
        NuANDNv = sum(ismember(Nu, Nv));
        NuMinusNv = length(setdiff(Nu, Nv));
        NvMinusNu = length(setdiff(Nv, Nu));
        
        Sfs(i, j)   = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
        Sfs(j, i) = Sfs(i, j);
        
        disp(i);
        disp(j);
    end
end

save Sfs.mat Sfs
