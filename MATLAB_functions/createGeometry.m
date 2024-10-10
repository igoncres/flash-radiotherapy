function [dl, R, coord] = createGeometry(geom, gridVessels)

%INFO: function to generate a voxel with a capillary distribution (circles)

%INPUT
%geom: [vector] parameters to create the voxel and vascular network
    %[voxel side length (um), vf, mean diameter (um), var diameter (um),...
    % min diameter (um), max diameter (um), intervessel distance (um)]
%gridVessels: 1 for a grid of vessels; otherwise random vessels

%OUTPUT
%dl: [matrix] geometry matrix (minimal regions)
%R: [vector] capillary radii
%coord: [matrix] spatial coordinates of capillary centers X(col 1) Y(col 2)

    %voxel size and vascular fraction
    V = geom(1);
    vf = geom(2);
    % mean and variance lognormal
    m = geom(3);
    v = geom(4);
    % max-min diameter of capillaries
    min_value = geom(5);
    max_value = geom(6);
    % min intervessel distance
    min_dist = geom(7);
    
	%create capillaries
    if gridVessels == 1
        [coord,R] = createVesselDistributionGrid(vf,V,m);
    else
        [coord,R] = createVesselDistribution(vf,V,m,v,max_value,min_value,min_dist);
    end

	%geometry matrix to decsg_fix
    [gd,sf,ns] = generateDescGeom(coord,V,R);

    w = size(coord,1);
    if w > 150
        fprintf('WARNING, number of capillaries may be too large, w= %i',w);
    end
 
    %generate dl matrix
    [dl] = decsg_fix(gd,sf,ns); 
    disp('dl matrix generated!')

    %patch I: if decsg_fix cannot generate dl, try again.
    while isempty(dl)
        disp('Empty dl, try again!')
        if gridVessels == 1
            [coord,R] = createVesselDistributionGrid(vf,V,m);
        else
            [coord,R] = createVesselDistribution(vf,V,m,v,max_value,min_value,min_dist);
        end
        [gd,sf,ns] = generateDescGeom(coord,V,R);
        [dl] = decsg_fix(gd,sf,ns);
    end
    
    %patch II: if decsg_fix misses any capillary, try again
    l_value=V;
    R_dl=[];
    for i=5:size(dl,2)
        if dl(10,i)~=l_value
            l_value=dl(10,i);
            R_dl=[R_dl l_value];
        end
    end
    k=0;
    
    while length(R)~=length(R_dl)
        k=k+1;
        if gridVessels == 1
            [coord,R] = createVesselDistributionGrid(vf,V,m);
        else
            [coord,R] = createVesselDistribution(vf,V,m,v,max_value,min_value,min_dist);
        end
        [gd,sf,ns] = generateDescGeom(coord,V,R);

        [dl] = decsg_fix(gd,sf,ns);
        while isempty(dl)
            if gridVessels == 1
                [coord,R] = createVesselDistributionGrid(vf,V,m);
            else
                [coord,R] = createVesselDistribution(vf,V,m,v,max_value,min_value,min_dist);
            end
            [gd,sf,ns] = generateDescGeom(coord,V,R);
            [dl] = decsg_fix(gd,sf,ns);
        end
        
        l_value=V;
        R_dl=[];
        for i=5:size(dl,2)
            if dl(10,i)~=l_value
                l_value=dl(10,i);
                R_dl=[R_dl l_value];
            end
        end
        
    end
end

%%

function [coord,R] = createVesselDistribution(vf,V,m,v,max_value,min_value,min_dist)

%INFO: function to generate vessels

%INPUT
%[vf, voxel side length (um), mean diameter (um), var diameter (um), ...
% max diameter (um), min diameter (um), intervessel distance (um)]

%OUTPUT
%coord: [matrix] spatial coordinates of capillary centers X(col 1) Y(col 2)
%R: [vector] capillary radii

    %Get the first radius
    R = getR(m,v);

    %Check if it fits the minimum and maximum diameter constraints
    i=1;
    while (R(i) < min_value/2 || R(i) > max_value/2)
        R = getR(m,v);
    end

    %Generate spatial coordinates, avoiding overlaps
    coord(1,1)=rand*(V-2*R(1)*1.001)+R(1)*1.001;
    coord(1,2)=rand*(V-2*R(1)*1.001)+R(1)*1.001;

    aa=0;
    bb=0;
    vf_min=(min_value/2)^2*pi/(V^2);
    vf_aux=R(1)*R(1)*pi/V^2;

    while vf_aux <= vf-vf_min || aa==1
        bb=0;
        i=i+1;

        R(i) = getR(m,v);

        while (R(i) < min_value/2 || R(i) > max_value/2)
            R(i) = getR(m,v);
        end

        coord(i,1)=rand*(V-2*R(i)*1.001)+R(i)*1.001;
        coord(i,2)=rand*(V-2*R(i)*1.001)+R(i)*1.001;

        distance1 = zeros(1,i-1);
        for j=1:i-1
            distance1(j)=norm(coord(i,:)-coord(j,:));
        end
        distance2=distance1-R(i)-R(1:i-1);

        if sum(distance2 < min_dist) > 0
            bb=1;
            i=i-1;    
        end
        clear distance1

        if bb==0
            vf_aux = sum(R.*R*pi)/V^2;
            if vf_aux > vf*1.1
                i=i-1;
                aa=1;
            else
                aa=0;
            end
        end

    end

end

%%

function [coord,R] = createVesselDistributionGrid(vf,V,m)

%INFO: function to generate the vessels

%INPUT
%[vf, voxel side length (um), diameter (um)]

%OUTPUT
%coord: [matrix] spatial coordinates of capillary centers X(col 1) Y(col 2)
%R: [vector] capillary radius

    vf_cap = pi*(m/2)^2/V^2;
    ncap = vf/vf_cap;
    
    ncap_row = round(sqrt(ncap));
    cap_dist = V/ncap_row;
    
    R(1:ncap_row^2) = m/2+rand(1,ncap_row^2)*0.1;
    
    coord = zeros(ncap_row^2,2);
    
    x=(cap_dist/2:cap_dist:V);
    y=(cap_dist/2:cap_dist:V);
    
    i=1;
    for j = 1:ncap_row
        for k = 1:ncap_row
            coord(i,1) = x(j);
            coord(i,2) = y(k);
            i=i+1;
        end
    end

end

%%

function [gd,sf,ns] = generateDescGeom(coord,V,R)

    N = size(coord,1);

    gd=zeros(10,N+1);
    gd(1,1)=3;
    gd(2,1)=4;
    gd(3,1)=0;
    gd(4,1)=V;
    gd(5,1)=V;
    gd(6,1)=0;
    gd(7,1)=0;
    gd(8,1)=0;
    gd(9,1)=V;
    gd(10,1)=V;  
   
    for i = 2:N+1    
        for j=1:4
            if j==1
                gd(j,i)=1;
            elseif j==2
                gd(j,i)=coord(i-1,1);
            elseif j==3
                gd(j,i)=coord(i-1,2);
            elseif j==4
%             gd(j,i)=R; %r=cte
                gd(j,i)=R(i-1);  %r=variable
            end
        end
    end
    
    sf='R1';
    for i=2:N+1
        cont=['-C' int2str(i-1)];
        sf = [sf cont];
    end

    ns = cell(N,1);
    for i = 1:N+1
        if i == 1
            ns(i,1) = cellstr(['R' int2str(i)]);
        else
            ns(i,1) = cellstr(['C' int2str(i-1)]);
        end
    end

    ns = char(ns);
    ns = ns';

end

%%

function R = getR(m,v)

    R = lognrnd(log(m),v)/2; %lognormal
%     R = exprnd(m)/2; %exponential

end