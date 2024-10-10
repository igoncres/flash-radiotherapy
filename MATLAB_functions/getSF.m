function [SF_conv, SF_flash] = getSF(u_flash, u, alpha, beta, alpha_beta, OER_a, OER_b, Km, D, model, T, dt)

%INFO: function to calculate the SFs with CONV-RT and FLASH-RT

%INPUT
%u_flash: [matrix] oxygenation at each node during FLASH-RT (mmHg)
%u: [vector] oxygenation at each node steady-state (mmHg)
%alpha, beta, alpha_beta: [real] LQ parameters and alpha/beta ratio
%OER_a, OER_b, Km: [real] OERs parameters
%D: [real] dose: flash(1)
%model: [struct] PDE model
%T: [real] time of dose delivery: flash(4)
%dt: [real] time step

%OUTPUT
%SF_conv: [real] surviving fraction with CONV-RT
%SF_flash: [real] surviving fraction with FLASH-RT

%Surviving fraction CONV-RT
[~, tri_areas, var_tri, ~, total_area]  = ponderation(u, model);
SF_conv = sum(calculateSF(var_tri, alpha, beta, D, Km, OER_a, OER_b).*tri_areas)/total_area;   

%Give BED as an input or calculate it
ntri = length(var_tri);
nt = size(u_flash,2);
BED(1:ntri,1:nt)=0;
for i=1:nt
    [~, ~, var_tri, ~, ~]  = ponderation(u_flash(:,i), model);
    BED(:,i) = getBED(var_tri, Km, OER_a, alpha_beta, D);
end

%Surviving fraction FLASH-RT
SF_flash = calculateFlashSF(BED, alpha, T, model, dt);

end

%%

function BED = getBED(u, K, OERmax, alpha_beta, D)

%INFO: Function to calculate the BED

%INPUT
%u: [vector] Oxygen level at the centroids of the mesh elements
%K: [real] Value at which the OER is a half of its maximum value
%OERmax: [real] OER value
%alpha/beta: [real] LQ parameters ratio
%D: [real] dose

%OUTPUT
%BED: [vector] BED at each node of the mesh

    OER = (OERmax.*u+K)./(u+K);
    BED = D*OER./OERmax.*(1+D*OER./(alpha_beta*OERmax));

end