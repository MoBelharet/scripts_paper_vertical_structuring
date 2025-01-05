function [ mod_, quartiles ] = estimate_statistic_metrics(V, Y)
% size(V,2) = 100
% length(Y) = 100 . Y est un vecteur

% normalization
V = V/nansum(V);

V_cum = cumsum(V,2);

[~,id_max_d] = max(V,[],2);


mod_ = Y(id_max_d);


q = [0.25 0.50 0.75];

for i=1:length(q)
    [~,id_q] = min(abs(V_cum - q(i)),[],2);
    quartiles(:,i) = Y(id_q);
    
end


end

