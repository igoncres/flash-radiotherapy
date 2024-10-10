function [SF_conv, SF_flash, TCP_conv, TCP_flash] = calculateTCP_het(geom, flash, dl, R, coord_cap, alpha_ox, alpha_beta, D, Km, OER_a, OER_b, dt, N0)

%INFO: function to calculate the TCP in heterogeneously oxygenated tumors

%INPUT
%geom, flash, dl, R, coord_cap: inputs to solve the oxygenation problem
%alpha_ox, alpha_beta: [real] LQ parameter and alpha/beta ratio
%D: [real] dose: flash(1)
%OER_a, OER_b, Km: [real] OERs parameters
%dt: [real] time step
%N0: [integer] number of cells

%OUTPUT
%SF_conv: [real] surviving fraction with CONV-RT
%SF_flash: [real] surviving fraction with FLASH-RT
%TCP_conv: [real] TCP with CONV-RT
%TCP_flash: [real] TCP with FLASH-RT

beta_ox=alpha_ox/alpha_beta;

%Oxygenation
[u0, u_flash, ~, ~, ~, model] = mainFLASH(geom, flash, dt, dl, R, coord_cap);

%SF
[SF_conv, SF_flash] = getSF(u_flash, u0, alpha_ox, beta_ox, alpha_beta, OER_a, OER_b, Km, D, model, flash(4), dt);

TCP_conv = exp(-N0*SF_conv);
TCP_flash = exp(-N0*SF_flash);

end