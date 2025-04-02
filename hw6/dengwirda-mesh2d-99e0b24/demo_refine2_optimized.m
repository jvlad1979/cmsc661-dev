function demo_refine2_optimized
% DEMO_REFINE2_OPTIMIZED - Mesh a semicircular domain with two circular holes
% using refined edge generation with reduced points to save memory.

    fprintf(1, [ ...
    ' Generating a mesh for a semicircular domain with two circular holes,', ...
    ' followed by refinement using refine2 (optimized version for memory usage).\n']);

    % Define outer semicircle with fewer points
    theta = linspace(0, pi, 15)';  % Reduced points for the semicircle
    outer_circle = [cos(theta), sin(theta)];

    % Define two circular holes with fewer points
    theta_hole = linspace(0, 2*pi, 10)';  % Reduced points for holes
    hole1 = 0.2 * [cos(theta_hole) - 0.25, sin(theta_hole) - 0.2];
    hole2 = 0.2 * [cos(theta_hole) + 0.25, sin(theta_hole) - 0.2];

    % Combine all nodes into one array
    node = [outer_circle; hole1; hole2];

    % Edge definitions
    n_outer = size(outer_circle,1);
    n_hole1 = size(hole1,1);
    n_hole2 = size(hole2,1);

    % Connect the outer semicircle
    edge_outer = [(1:n_outer-1)', (2:n_outer)'; n_outer, 1];

    % Connect hole 1
    start_h1 = n_outer + 1;
    edge_hole1 = [(start_h1:start_h1+n_hole1-2)', (start_h1+1:start_h1+n_hole1-1)';
                  start_h1+n_hole1-1, start_h1];

    % Connect hole 2
    start_h2 = start_h1 + n_hole1;
    edge_hole2 = [(start_h2:start_h2+n_hole2-2)', (start_h2+1:start_h2+n_hole2-1)';
                  start_h2+n_hole2-1, start_h2];

    % Combine all edges
    edge = [edge_outer; edge_hole1; edge_hole2];

    % Perform mesh refinement with a larger target edge length
    hfun = 0.2;  % Larger edge length to reduce memory usage

    % Call refine2 with the optimized edge and node set
    [vert, etri, tria, tnum] = refine2(node, edge, [], [], hfun);

    % Plot the refined mesh
    figure;
    patch('faces', tria(:,1:3), 'vertices', vert, 'facecolor', 'w', 'edgecolor', [.2,.2,.2]);
    hold on;
    plot(node(:,1), node(:,2), 'ro', 'markerfacecolor', 'r');
    title('Refined Mesh for Semicircle with Two Holes');
    axis equal;

end
