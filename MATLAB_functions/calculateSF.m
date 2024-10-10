function [SF, alpha, beta] = calculateSF(u, alpha, beta, D, K, OER_a, OER_b)

%INFO: function to evaluate the surviving fraction in CONV-RT with OERs

%INPUT
%u: [vector] oxygenation at the centroids of the mesh elements
%alpha: [real] alpha parameter LQ
%beta: [real] beta parameter LQ
%D: [real] dose 
%K, OER_a, OER_b: [real] standard parameters for the calculation of OERs 

%OUTPUT
%SF: [vector] surviving fraction at each node

    [OER_a, OER_b] = getOERs(u, K, OER_a, OER_b);
    alpha = alpha./OER_a;
    beta = beta./(OER_b.^2);
    SF = exp(-alpha.*D - beta.*(D.^2));

end

%%

function [OER_a, OER_b] = getOERs(u, K, OER_ma, OER_mb)

%INFO: function to calculate the OERs (Sovik et al., 2006)

    OER_a = OER_ma.*(u+K)./((u.*OER_ma)+K);
    OER_b = OER_mb.*(u+K)./((u.*OER_mb)+K);

end