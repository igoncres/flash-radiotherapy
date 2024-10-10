function [param_best, cost_best] = annealing_optimization(param_0, T0, dT, nit_max, nit_T, u, u_flash, D, model, T, dt, tspan, RTtime, exp_times, exp_vol, exp_err)

%INFO: Simulated annealing to fit model to experimental data

%INPUT
%param_0: [vector] initial parameters
%T_0: [real] initial temperature
%dT: [real] temperature delta
%nit_max: [integer] number of steps in which temperature decreases
%nit_T: [integer] number of evaluations for each temperature
%u: [vector] oxygen at each node of the mesh
%u_flash: [vector] oxygen at each node of the mesh during FLASH-RT
%D: [real] dose
%model: [struct] PDE model
%T: [real] time of irradiation
%dt: [real] time step in the oxygenation problem
%tspan: [vector] initial and final time of volume simulation
%RTtime: [vector] times of radiation delivery
%exp_times: [struct] times of experimental measurements
%exp_vol: [struct] experimental volumes
%exp_err: [struct] experimental uncertainties

%OUTPUT
%param_best: [vector] optimal parameters
%cost_best: [real] optimal cost

% Parameters and cost initialization
    param_op = param_0;
    param_best = param_0;
    cost_op = least_squares(param_0, u,  u_flash, D, model, T, dt, tspan, RTtime, exp_times, exp_vol, exp_err);
    cost_best = cost_op;

    Temp = T0; % Initial temperature

    n = nit_max; % Number of steps in which temperature decreases
    m = nit_T;   % Number of evaluations for each temperature

    for i = 1:n
        for j = 1:m
            j
            % Calculate new parameters
            param_new = neighbor(param_op, Temp, T0);
            
            % Calculate the cost associated to the new parameters
            cost_new = least_squares(param_new, u,  BED, D, model, T, dt, tspan, RTtime, exp_times, exp_vol, exp_err);
            
            % Calculate the survival for the new parameters
            surv = survival(cost_new, cost_op, Temp);
            if surv == 1
                if cost_new < cost_best
                    cost_best = cost_new;
                    param_best = param_new;
                end
           
                param_op = param_new;
                cost_op = cost_new;
            end
        end %for j
    
        % cooling
        Temp = dT * Temp;

    end %for i

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%     Function to obtain neighbor parameters    %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = neighbor(param, T, T0)

    var = (2 * rand - 1.0) * sqrt(T / T0);
    m = length(param);    
    id = ceil(rand * m);

    param(id) = max(0, param(id) * (1 + var));

    %Constraints
    param(1) = min(log(2), max(log(2)/20, param(1)));           % lambda
    param(2) = min(1, max(0, param(2)));                        % phi
%     param(3) = min(1e5, max(1500, param(3)));                 % K Diffenderfer (a,b) & Zhu (d)
    param(3) = min(1e5, max(300, param(3)));                    % K Zhu (d)
    param(4) = min(1, max(0.01, param(4)));                     % alpha
    param(5) = min(21.5+16, max(21.5-16, param(5)));            % V0 Zhu (d)
%     param(6) = min(768, max(202, param(6)));                  % V0 Zhu (c) Control
%     param(7) = min(255, max(65, param(7)));                   % V0 Zhu (c) flash
%     param(8) = min(507, max(217, param(8)));                  % V0 Zhu (c) conv

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%     Function to asign a suvival probability to a solution     %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function surv = survival(cost_new, cost, T)

    surv = 0;

    if cost_new < cost
        surv = 1;
    elseif rand < exp(- (cost_new - cost) / T)
        surv = 1;
    end
    
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        %
%     Objective function computation     %
%                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = least_squares(param_new, u,  u_flash, D, model, T, dt, tspan, RTtime, exp_times, exp_vol, exp_err)

    cost = fittingFunction(param_new, u, u_flash, D, model, T, dt, tspan, RTtime, exp_times, exp_vol, exp_err);

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %
%     Fitting function    %
%                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = fittingFunction(par, u, u_flash, D, model, T, dt, tspan, RTtime, exp_times, exp_vol, exp_err)

lambda = par(1);
phi = par(2);
K = par(3);
alpha = par(4);
V0 = par(5);
% V0 =par(6);
% V0F = par(7);
% V0C = par(8);
beta = alpha/10;
OER_a = 2.5;
OER_b = 2.5;
Km = 3.28;
dt_vol = 0.1;

[SF_conv, SF_flash] = getSF(u_flash, u, alpha, beta, alpha_beta, OER_a, OER_b, Km, D, model, T, dt);

% V0=485;
% V0C=362;
% V0F=210;

[~,~,~,V_notreat] = getVolume(V0, 1, lambda, 0, K, tspan, RTtime, dt_vol);
[~,~,~,V_conv] = getVolume(V0, SF_conv, lambda, phi, K, tspan, RTtime, dt_vol);
[~,~,~,V_flash] = getVolume(V0, SF_flash, lambda, phi, K, tspan, RTtime, dt_vol);

id1 = cell2mat(exp_times(1,:))./0.1+1;
fnotreat = norm((V_notreat(id1)-cell2mat(exp_vol(1,:)))./cell2mat(exp_err(1,:)));
id2 = cell2mat(exp_times(2,:))./0.1+1;
fconv = norm((V_conv(id2)-cell2mat(exp_vol(2,:)))./cell2mat(exp_err(2,:)));
id3 = cell2mat(exp_times(3,:))./0.1+1;
fflash = norm((V_flash(id3)-cell2mat(exp_vol(3,:)))./cell2mat(exp_err(3,:)));

f = fnotreat+fconv+fflash;

end