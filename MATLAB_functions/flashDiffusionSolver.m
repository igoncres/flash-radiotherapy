function [u, sol] = flashDiffusionSolver(model, T, G0, k_ROD, D, tlist)

%INFO: function to solve the O2 diffusion problem with UHDR

%INPUT
%model: [struct] PDE model
%T: time of dose delivery
%G0: intrinsic rate of radiolytic oxygen depletion (ROD)
%k_ROD: oxygen partial pressure at which ROD is half its maximum
%D: RT dose

%OUTPUT
%u: [vector] pO2 distribution
%sol: [struct] solution info
    
    % Specify coefficients, diffusion coefficient as c
    acoeff = @(location,state) acoeffunction(location, state, D, G0, T, k_ROD);
    D = 2000; %Oxygen diffusion coefficient
    specifyCoefficients(model, 'Face', 1, 'm', 0, 'd', 1, 'c', D, 'a', acoeff, 'f', 0);
    
    % Solve
    sol = solvepde(model,tlist);
    
    u = sol.NodalSolution;
    u(u<0) = 0;
    
end

%%

function acoeff = acoeffunction(~, state, D, G0, T, k_ROD)

    state_u = state.u;
    acoeff = 15*(1./(state_u + 2.5)) + D/T*G0 * (1./(state_u + k_ROD));

end