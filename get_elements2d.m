function [e] = get_elements2d(mesh)
    % function GET_ELEMENTS2D - 2D elements of a 2D or 3D mesh.
    
    if isa(mesh, 'mesh2d')
        e = mesh.elements;
    else
        e = mesh.elements2d;
    end
end
