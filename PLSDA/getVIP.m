function vip = getVIP(p, ncomp, W, pctVar)
% get VIP values
vip = zeros(p,1);
% sumW = sum(W);
sumW = sqrt(sum(W.^2));
W = W./repmat(sumW, p, 1);
for i = 1:p
    a = sum(W(i,:).^2 .* pctVar);
    b = sum(pctVar);
    tmp = sqrt(p*a/b);
    vip(i) = tmp;
%     t = sqrt(p/sum(tmp.^2));
%     vip(i) = tmp * t;
    
    
%     vip(i) = sqrt(p*a);
end


end

