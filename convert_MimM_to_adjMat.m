load Mim5NN;

MimM_ADJ = zeros(5080);

for i = 1 : 5080
    MimM_ADJ(i, :) = MimM(i, :) >= 0.2;
end

save MimAdj.mat MimM_ADJ