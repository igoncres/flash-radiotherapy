function [u, u_flash, R, coord_cap, dl, model, time, pO2, sol] = mainFLASH(geom, flash, dt, varargin)

%INFO: function to obtain pO2 distribution before, during and after FLASH-RT

%INPUT
%geom: [vector] parameters that describe the vascular geometry
    %[V, vf, m, v, min_value, max_value, int_dist] 
    %V: voxel side (um)
    %vf: vascular fraction (0,1)
    %[m, v]: lognormal parameters (um)
    %[min_value max_value]: min/max size of generated capillaries (um)
    %int dist: minimum intervessel distance (um)
%flash: [vector] parameters that describe FLASH-RT effect
    %[D G0 k_ROD T]
    %D: radiation dose
    %G0: intrinsic rate of radiolytic oxygen depletion (ROD)
    %k_ROD: oxygen partial pressure at which ROD is half its maximum
    %T: time of dose delivery = D/dose_rate
%varargin: variable number of input arguments (geometry can be externally
	%provided or generated within the function)
    %dl, R, coord_cap, model, u0

%OUTPUT
%u: [vector] steady-state pO2 distribution (initial)
%u_flash: [matrix: mesh nodes x time nodes] pO2 distribution during and/or after FLASH-RT
%R: [vector] radii of capillaries
%coord_cap: [vector] spatial coordinates of capillaries
%dl: [matrix] geometry matrix (minimal regions)
%model: [struct] PDE model
%time: [vector] time associated to the variable pO2
%pO2: [vector] mean pO2 before and after FLASH-RT
%sol: [struct] steady-state solution
        
%RT Parameters
    D = flash(1);
    G0 = flash(2);
    k_ROD = flash(3);
    T = flash(4);

%Geometry and mesh
    nVarargs = length(varargin);
    switch nVarargs
        case 3
        dl = varargin{1};
        R = varargin{2};
        coord_cap = varargin{3};
        model = [];
        u0 = [];
        case 4
        dl = varargin{1};
        R = varargin{2};
        coord_cap = varargin{3};
        model = varargin{4};
        u0 = [];
        case 5
        dl = varargin{1};
        R = varargin{2};
        coord_cap = varargin{3};
        model = varargin{4};
        u0 = varargin{5};
        otherwise
        dl = []; R = []; coord_cap = []; model = []; u0 = [];
    end

    if isempty(dl)
        [dl, R, coord_cap] = createGeometry(geom, 0);
    end

    if isempty(model)
        %Create PDE model
        model = createpde(1);
        %Assign geometry to model
        geometryFromEdges(model, dl);
        disp('Geometry!')
        %Create mesh
        generateMesh(model,'Hmin', 1, 'Hmax', 10, 'GeometricOrder', 'linear')
        disp('Mesh!')
    end

%Model
    if isempty(u0)
        % Set initial condition
%         delete(model.InitialConditions)
        setInitialConditions(model, 0);
    else
        setInitialConditions(model, u0);
    end
    
    % Identify edges with capillaries
    id_cap_edge = assignEdges(dl);

    % Assign boundary conditions
%     delete(model.BoundaryConditions)
    applyBoundaryCondition(model, 'neumann', 'Edge', 1:4, 'g', 0, 'q', 0);
    for i = 1:length(R)
        applyBoundaryCondition(model, 'dirichlet', 'Edge', id_cap_edge{i}, 'u', 40);
    end


%Solve O2 diffusion equation (steady-state solution)
    [u, sol] = pO2DiffusionSolver(model);
    pO2(1) = ponderation(u, model);
    time(1) = 0;

%     figure();
%     pdeplot(model,'XYData',u);
    
%Solve the oxygenation problem with ROD
    delete(model.InitialConditions)
    setInitialConditions(model, sol); %steady-state solution as initial condition
    tlist = (0:dt:T);
    u_flash = flashDiffusionSolver(model, T, G0, k_ROD, D, tlist);
    pO2(2) = ponderation(u_flash(:,end), model);
    time(2) = T;
    
end

%%

function id_cap_edge = assignEdges(dl)

%INFO: function to assign mesh edges to each vessel

%INPUT
%dl: [matrix] mesh matrix info

%OUTPUT
%id_cap_edge: [cell] edge ids for each vessel

    id_cap_edge{1} = 5; %First capillary edge
    nedges = length(dl);
    j = 1;
    for i = 5:nedges-1
        if abs(dl(10,i+1)-dl(10,i)) < 1e-10
            id_cap_edge{j} = [id_cap_edge{j} i+1];
        else
            j = j+1;
            id_cap_edge{j} = i+1;
        end
    end

end