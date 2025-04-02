function [vert, tria, edges] = distmesh_refine2(fd, fh, hmin, box, nodes, edges)
% DISTMESH_REFINE2 - A wrapper for refine2 to behave like distmesh
%   fd      : distance function
%   fh      : mesh size function
%   hmin    : minimum edge length (for mesh refinement)
%   box     : bounding box [xmin, xmax; ymin, ymax]
%   nodes   : initial nodes (optional)
%   edges   : initial edges (optional)

    % If no nodes or edges are provided, generate initial ones
    if nargin < 5 || isempty(nodes)
        nodes = generate_nodes(fd, box, fh);
    end
    if nargin < 6 || isempty(edges)
        edges = generate_edges(nodes, fd, box);
    end

    % Create the initial mesh (nodes, edges, triangles)
    [vert, etri, tria, tnum] = refine2(nodes, edges);

    % Refine the mesh based on the desired edge length
    [vert, etri, tria, tnum] = refine2(vert, etri, tria, tnum, hmin);

    % Return the refined mesh: vertices, triangles, and edges
    edges = generate_edges(vert, fd, box);
end

% Generate the nodes by sampling the distance function
function nodes = generate_nodes(fd, box, fh)
    [X, Y] = meshgrid(box(1,1):0.02:box(1,2), box(2,1):0.02:box(2,2)); % Grid sampling
    nodes = [X(:), Y(:)];

    % Filter nodes by the distance function (inside the domain)
    valid_points = fd(nodes) < 0;
    nodes = nodes(valid_points, :);
end

% Generate edges by connecting nearby nodes
function edges = generate_edges(nodes, fd, box)
    edges = [];
    for i = 1:size(nodes, 1)
        for j = i+1:size(nodes, 1)
            dist = sqrt((nodes(i,1)-nodes(j,1))^2 + (nodes(i,2)-nodes(j,2))^2);
            if dist < 0.03  % A threshold to connect nearby nodes
                edges = [edges; i, j];
            end
        end
    end
end

