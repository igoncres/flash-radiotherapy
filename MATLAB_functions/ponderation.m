function [var_pond, tri_areas, var_tri, center_tri, total_area]  = ponderation(v, model)

%INFO: function to ponderate any variable per node taking into account the 
    %different weight of each node (triangle area)
    
%INPUT:
%v: [vector] value of the variable at each node
%model: [struct] PDE model

%OUTPUT:
%var_pond: [real] ponderated value
%tri_areas: [vector] triangle areas
%var_tri: [vector] average value per triangle
%center_tri: [matrix] triangle center position

    modelMesh = model.Mesh;
    meshElements = modelMesh.Elements;
    meshNodes = modelMesh.Nodes;

    %number of triangles
    ntri = length(meshElements);

    var_tri = zeros(1, ntri);
    center_tri = zeros(2, ntri);

    for k = 1:ntri
        tri_points_id = meshElements(:,k); %Points id
        tri_points = meshNodes(:,tri_points_id); %Points coordinates
        center_tri(:,k) = mean(tri_points,2); %Triangle centroids
        %Non-linear elements
%         var_tri(k) = interpolateSolution(sol,center_tri(:,k));
        %Linear elements
        var_tri(k) = mean(v(tri_points_id));
    end
    
    [total_area, tri_areas] = area(modelMesh);
    var_pond = sum(tri_areas.*var_tri)/total_area;
    
end