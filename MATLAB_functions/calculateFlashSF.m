function SF = calculateFlashSF(BED, alpha, T, model, dt)

%INFO: evaluate the surviving fraction in FLASH-RT (Taylor et al., 2022)

%INPUT
%BED: [vector] BED at the centroids of the mesh elements
%alpha: [real] alpha parameter LQ
%T: [real] final time
%model: [struct] PDE model
%dt: [real] time step

%OUTPUT
%SF: [real] surviving fraction

    nt = size(BED,2);
    int_p = zeros(1,nt);
    [areaMesh, tri_areas] = area(model.Mesh);

    for i=1:nt
        aux = exp(-alpha*BED(:,i));
        int_p(i) = sum(aux.*tri_areas'./areaMesh);
    end

    int_t = sum(int_p*dt);
    SF = 1/T*int_t;
    
end