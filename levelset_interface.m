function [contour, element_indices] = levelset_interface(md, levelset)
    % Function LEVELSET_INTERFACE - zero isocontour of the levelset
    %
    % The contour where the levelset field is zero as an unordered list of vertex pairs.
    % Each vertex pair is the points where the contour intersects an element of the mesh, forming a line segment of the contour.
    % Intersections are computed by linear interpolation of the levelset field along the edges of an element.
    % Returns an N x 2 x 2 matrix where N is the number of elements that are intersected by the contour.
    % contour(i, 1, 1) and contour(i, 1, 2) are the x and y coordinates of the first vertex of the i-th contour segment.
    % contour(i, 2, 1) and contour(i, 2, 2) are the x and y coordinates of the second vertex of the i-th contour segment.
    %
    % md: ISSM model, 2d or 3d, (extruded) triangle mesh.
    % levelset: levelset field to find the zero isocontour for.
    %

    % which elements are intersected by the zero isocontour
    mesh = md.mesh;
    elements = get_elements2d(mesh);
    n_in = sum(levelset(elements) <= 0, 2);
    is_interface = n_in > 0 & n_in < 3;
    interface_elements = elements(is_interface, :);
    
    % compute points where each element is intersected
    num_elements = size(interface_elements, 1);
    contour = zeros(num_elements, 2, 2);
    if nargout > 1
        element_indices = find(is_interface);
    end
    for i = 1:num_elements
        current_element = interface_elements(i, :);

        %find the edges that are intersected by the interface
        vert_in = current_element(levelset(current_element) <= 0);
        vert_out = current_element(levelset(current_element) > 0);
        [I, O] = meshgrid(vert_in, vert_out);
        interface_edges = [I(:), O(:)];
        
        %pt 1
        edge1 = interface_edges(1, :);
        v1 = abs(levelset(edge1(1)));
        v2 = abs(levelset(edge1(2)));
        w = v1 / (v1 + v2);
        pt1 = w * [mesh.x(edge1(2)), mesh.y(edge1(2))] + (1-w) * [mesh.x(edge1(1)), mesh.y(edge1(1))];

        %pt 1
        edge2 = interface_edges(2, :);
        v1 = abs(levelset(edge2(1)));
        v2 = abs(levelset(edge2(2)));
        w = v1 / (v1 + v2);
        pt2 = w * [mesh.x(edge2(2)), mesh.y(edge2(2))] + (1-w) * [mesh.x(edge2(1)), mesh.y(edge2(1))];
        
        contour(i, 1, :) = pt1;
        contour(i, 2, :) = pt2;
    end
end
