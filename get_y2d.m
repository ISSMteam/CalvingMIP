function [e] = get_y2d(mesh)
    % function GET_Y2D - 2D y coordinates of a 2d or 3d mesh
    
    if isa(mesh, 'mesh2d')
        e = mesh.y;
    else
        e = mesh.y2d;
    end
end
