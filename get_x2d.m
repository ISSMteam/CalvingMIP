function [e] = get_x2d(mesh)
    % function GET_X2D - 2D x coordinates of a 2D or 3D mesh.
    if isa(mesh, 'mesh2d')
        e = mesh.x;
    else
        e = mesh.x2d;
    end
end
