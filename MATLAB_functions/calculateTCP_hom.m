function [SF_conv, SF_flash, TCP_conv, TCP_flash] = calculateTCP_hom(u0, alpha_ox, alpha_beta, D, dose_rate, Km, OER_a, OER_b, G0, k_ROD, dt, N0)

%INFO: function to calculate the TCP in homogeneously oxygenated tumors

%INPUT
%u0: [real] oxygenation at each node during FLASH-RT (mmHg)
%alpha_ox, alpha_beta: [real] LQ parameter and alpha/beta ratio
%D: [real] dose: flash(1)
%OER_a, OER_b, Km: [real] OERs parameters
%G0: [real] radiolytic consumption rate
%k_ROD: [real] oxygen level at which the radiolytic consumption is half of its maximum
%dt: [real] time step
%N0: [integer] number of cells

%OUTPUT
%SF_conv: [real] surviving fraction with CONV-RT
%SF_flash: [real] surviving fraction with FLASH-RT
%TCP_conv: [real] TCP with CONV-RT
%TCP_flash: [real] TCP with FLASH-RT

beta_ox=alpha_ox/alpha_beta;

%CONV-RT
SF_conv = calculateSF(u0, alpha_ox, beta_ox, D, Km, OER_a, OER_b);

%FLASH-RT
T = D/dose_rate;
tspan = (0: dt: T);

[t,u] = ode45(@(t,u) -D/T*G0*u/(k_ROD+u), tspan, u0);

SF = t;
alpha = t;
beta = t;
for i=1:length(t)
    [SF(i), alpha(i), beta(i)] = calculateSF(u(i), alpha_ox, beta_ox, D, Km, OER_a, OER_b);
end
SF_flash = mean(SF);

TCP_conv = exp(-N0*SF_conv);
TCP_flash = exp(-N0*SF_flash);

end