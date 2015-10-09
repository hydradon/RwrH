function [B1,B2] = getBridgeM(bridgeM)
idx3 = find(sum(bridgeM) > 0);    
%     idx4 = setdiff([1:Ng],idx3);
    B1 = bridgeM; B1(:) = 0;    
    for ii = 1 : length(idx3)
        B1(:,idx3(ii)) = bridgeM(:,idx3(ii))/sum(bridgeM(:,idx3(ii)));        
    end
%     B1(:,idx4) = 1/Nd;

    B = bridgeM';
    idx5 = find(sum(B) > 0);
%     idx6 = setdiff([1:Nd],idx5);
    B2 = B; B2(:) = 0;   
    for ii = 1 : length(idx5)
        B2(:,idx5(ii)) = B(:,idx5(ii))/sum(B(:,idx5(ii))); 
    end