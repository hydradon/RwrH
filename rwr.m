function [PT,k] = rwr(W,P0,gamma)

    PT = P0;
    k = 0;
    delta = 1;
    while  (delta > 1e-10)
        PT1 = (1-gamma)*W*PT+gamma*P0;
        delta = sum(abs(PT1 - PT));
%         disp(delta);
        PT = PT1;
        k = k + 1;
    end