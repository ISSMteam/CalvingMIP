function [sum_area] = levelset_area(md, levelset, poly)
    %function LEVELSET_AREA - area of the region where the levelset is negative
    %
    % Sum of the area of elements of the mesh that are inside the region where the levelset is negative.
    % Partial area is computed for elements that are partially in the region by linear interpolation of the levelset field.
    % If a polygon is passed as the optional third argument, only the area inside that polygon is computed.
    % Partial area is computed for elements that are partially within the polygon by intersecting the element with the polygon. 
    % Assumes a mesh of (extruded) triangles.
    % 
    % md: ISSM model
    % levelset: levelset field that defines the area
    % poly: optional, polygon ([x1 y1; x2, y2; ...]) that defines the polygon
    %
    % e.g. ice_area = levelset_area(md, md.mask.ice_levelset)
    %      ice_area_northeast = levelset_area(md, md.mask.ice_levelset, [0 0; 1e6 0; 1e6 1e6; 0, 1e6])

    if nargin == 2
        % no polygon given, use polygon that encloses the whole model
        [minx, maxx] = bounds(md.mesh.x);
        [miny, maxy] = bounds(md.mesh.y);
        lx = maxx - minx;
        ly = maxy - miny;
        minx = minx - lx;
        maxx = maxx + lx;
        miny = miny - ly;
        maxy = maxy + ly;
        poly = [minx, miny; maxx, miny; maxx, maxy; minx, maxy];
    end

    [front_vertices, front_element_indices] = levelset_interface(md, levelset);
    el2d = get_elements2d(md.mesh);
    sum_area = 0;

    for i = 1:size(el2d, 1)
        % area of each element
        vertex_is_inside = levelset(el2d(i, :)) <= 0;
        if sum(vertex_is_inside) == 0
            continue; % no vertices in the region of negative levelset, ignore element
        end

        element_front_vertex_idx = find(front_element_indices == i);

        % find the polygon that describes the part of the element that is inside the negative levelset 
        if isempty(element_front_vertex_idx)
            % no front vertices -> whole element inside the front
            p = [md.mesh.x(el2d(i, :)), md.mesh.y(el2d(i, :))];
        else
            % one or two points inside the front ... 
            p = [md.mesh.x(el2d(i, find(vertex_is_inside))), md.mesh.y(el2d(i, find(vertex_is_inside)))];
            % ... plus the front segment form a triangle or quadrilateral.
            p = [p; flipud(squeeze(front_vertices(element_front_vertex_idx, :, :)))];
            % for quadrilateral, the order of the front vertices needs to be flipped, otherwise the polygon self-intersects.
            % (depends on the order of vertices delivered by the `levelset_interface` function)
            % for triangle, order doesn't matter.
        end

        %intersect with polygon
        vertex_is_inside = inpolygon(p(:, 1), p(:, 2), poly(:, 1), poly(:, 2));
        if ~any(vertex_is_inside)
            % all vertices outside the polygon, ignore element
            continue;
        elseif all(vertex_is_inside)
            % all vertices inside the polygon
            % cheap path, area of a polygon is quick
            sum_area = sum_area + polyarea(p(:, 1), p(:, 2)); 
        else
            % some vertices inside the polygon, use polygon intersection
            % expensive path, matlab functions `polyshape` and `intersect` are slow, but hopefully only a few elements.
            % probably could be replaced with a faster implementations.
            sum_area = sum_area + area(intersect(polyshape(p), polyshape(poly)));
        end
    end
end
