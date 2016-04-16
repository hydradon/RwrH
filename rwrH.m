% function [p,d,step]=rwrH(subPPI,subSfsGG,subMim,B1,B2,gamma,lamda,eta,alpha,d0,p0)
function [p,d,step]=rwrH(subPPI,subMim,B1,B2,gamma,lamda,eta,d0,p0)

Nd = size(subMim,1); %5080
Ng = size(subPPI,1); %8919
% ggmat = alpha * subPPI + (1-alpha) * subSfsGG;
% nW = [(1-lamda)*ggmat',lamda*B1'; lamda*B2',(1-lamda)*subMim']';
nW = [(1-lamda)*subPPI',lamda*B1'; lamda*B2',(1-lamda)*subMim']';
for i = 1 : size(nW,1)
    nW(:,i) = nW(:,i)/sum(nW(:,i));
end

P0 = [(1-eta)*p0;eta*d0];
% P0 = [p0;d0];
P0 = P0/sum(P0);
[PT,k] = rwr(nW,P0,gamma);
% disp(max(abs(PT-P0)))
p = PT(1:Ng);
d = PT(Ng+1:end);
step = k;