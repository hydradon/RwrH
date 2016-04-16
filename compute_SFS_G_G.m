load PPIM;

Ng = 8919;

Sfs_G_G = zeros(Ng);

for i = 1 : Ng
   disp(['i is' num2str(i)]);
   first_neighbors = find(PPIM(i, :));
   
   Nu = first_neighbors;
   if (PPIM(i, i) == 0)
       Nu(end + 1) = i;
   end
   
   %calculate Sfs distance between i and its first neighbors
   for j = 1 : length(first_neighbors)
       if first_neighbors(j) ~= i
           Nv = find(PPIM(first_neighbors(j), :));
           if (PPIM(first_neighbors(j), first_neighbors(j)) == 0)
               Nv(end + 1) = first_neighbors(j);
           end
           
           NuANDNv = sum(ismember(Nu, Nv));
           NuMinusNv = length(setdiff(Nu, Nv));
           NvMinusNu = length(setdiff(Nv, Nu));
           
           Sfs_G_G(i, first_neighbors(j))   = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
           Sfs_G_G(first_neighbors(j), i) = Sfs_G_G(i, first_neighbors(j));
       end
   end
   
   %calculate Sfs distance between i and its second neigbhors
   for j = 1 : length(first_neighbors)
      second_neighbors = find(PPIM(first_neighbors(j), :));

      for z = 1 : length(second_neighbors)
         Nv = find(PPIM(second_neighbors(z), :));
         if (PPIM(second_neighbors(z), second_neighbors(z)) == 0)
             Nv(end + 1) = second_neighbors(z);
         end

         NuANDNv = sum(ismember(Nu, Nv));
         NuMinusNv = length(setdiff(Nu, Nv));
         NvMinusNu = length(setdiff(Nv, Nu));

         if Sfs_G_G(i, second_neighbors(z)) ~= 0
             Sfs_G_G(i, second_neighbors(z))   = ((2*NuANDNv)/(NuMinusNv + 2*NuANDNv)) * ((2*NuANDNv)/(NvMinusNu + 2*NuANDNv));
             Sfs_G_G(second_neighbors(z), i) = Sfs_G_G(i, second_neighbors(z));
         end
      end
   end
end

save Sfs_G_G_2.mat Sfs_G_G