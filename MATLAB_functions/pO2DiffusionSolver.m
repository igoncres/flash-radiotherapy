function [u, sol] = pO2DiffusionSolver(model)

%INFO: function to solve the O2 diffusion problem

%INPUT
%model: [struct] PDE model

%OUTPUT
%u: [vector] solution
%sol: [struct] solution info
    
    % Specify coefficients, diffusion coefficient as c
    D = 2000; %Oxygen diffusion coefficient
    specifyCoefficients(model, 'Face', 1, 'm', 0, 'd', 0, 'c', D, 'a', @(location,state) acoeffunction(location, state), 'f', 0);
    
    % Solve
    sol = solvepde(model);
    
    u = sol.NodalSolution;
    u(u<0) = 0;
    
end

%%

function acoeff = acoeffunction(~, state)

    acoeff = 15 * (1 ./ (state.u + 2.5));
    
end