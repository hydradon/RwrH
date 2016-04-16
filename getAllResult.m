% results = zeros(9, 9, 9); %lamda x gamma x beta
T = table;
filename = 'testResult4.xlsx';
lower = 0.1;
upper = 0.9;
step = 0.1;
for lamda = 0.1:0.1:0.5
%     lamda = 0.9;
%     eta = 0.5;
%     gamma = 0.5;
%     disp(['lamda is ' num2str(lamda)]);
    for gamma = 0.1:0.1:0.9
%         disp(['gamma is ' num2str(gamma)]);
        for eta = 0.1:0.1:0.9
%             disp(['eta is ' num2str(eta)]);
            disp(['lamda is ' num2str(lamda) ' gamma is ' num2str(gamma) ' eta is ' num2str(eta)]);
            [LOO1, LOO2] = LOO(lamda, gamma, eta);
            disp(['LOO1 is ' num2str(LOO1) ' LOO2 is ' num2str(LOO2)]);
%             Tnew = struct('lamda',lamda,'gamma',gamma,'eta',eta,'LOO1',0,'LOO2',0);
            Tnew = struct('lamda',lamda,'gamma',gamma,'eta',eta,'LOO1',LOO1,'LOO2',LOO2);
            T = [T; struct2table(Tnew)];
        end
    end
end
writetable(T, filename);
disp(T)