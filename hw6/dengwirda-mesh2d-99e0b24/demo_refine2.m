function demo_refine2
% DEMO_REFINE2 - Mesh different geometries using refine2

    fprintf(1, 'Generating meshes for different geometries using refine2.\n');

    % Run different test cases
    generate_mesh_L();
    generate_mesh_pentagons();
    generate_mesh_semicircle();
end

function generate_mesh_L()
    fprintf(1, 'Generating L-shaped region mesh.\n');

    node = [
        0, 0; 1, 0; 1, 0.4; 0.4, 0.4; 0.4, 1; 0, 1  % L shape vertices
    ];

    edge = [
        1, 2; 2, 3; 3, 4; 4, 5; 5, 6; 6, 1 % L shape edges
    ];

    generate_and_plot_mesh(node, edge, 'L-shaped Region');
end

function generate_mesh_pentagons()
    fprintf(1, 'Generating pentagon inside a pentagon mesh.\n');

    outer_pv = 0.5 + [cos((0:4)'*2*pi/5 + pi/2), sin((0:4)'*2*pi/5 + pi/2)];
    inner_pv = 0.5 * ([0.5+cos((0:4)'*2*pi/5), 0.5+sin((0:4)'*2*pi/5)] - 0.5) + 0.5;

    node = [outer_pv; inner_pv];

    edge = [
        1, 2; 2, 3; 3, 4; 4, 5; 5, 1; % Outer pentagon edges
        6, 7; 7, 8; 8, 9; 9, 10; 10, 6 % Inner pentagon edges
    ];

    generate_and_plot_mesh(node, edge, 'Pentagon Inside Pentagon');
end

function generate_mesh_semicircle()
    fprintf(1, 'Generating semicircle with two holes mesh.\n');

    % Define outer semicircle
    theta = linspace(0, pi, 8)';  % More points for a smoother semicircle
    outer_circle = [cos(theta), sin(theta)];

    % Define two circular holes
    theta_hole = linspace(0, 2*pi, 8)';  % Full circles for holes
    hole1 = 0.2 * [cos(theta_hole) - 0.25, sin(theta_hole) - 0.2];
    hole2 = 0.2 * [cos(theta_hole) + 0.25, sin(theta_hole) - 0.2];

    % Combine all nodes
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

    generate_and_plot_mesh(node, edge, 'Semicircle with Holes');
end

function generate_and_plot_mesh(node, edge, title_str)
    [vert, etri, tria, tnum] = refine2(node, edge);

    figure;
    patch('faces', tria(:,1:3), 'vertices', vert, 'facecolor', 'w', 'edgecolor', [.2,.2,.2]);
    hold on;
    plot(node(:,1), node(:,2), 'ro', 'markerfacecolor', 'r');
    title(['Initial Mesh: ', title_str]);
    axis equal;

##    hfun = 0.1;
##    [vert, etri, tria, tnum] = refine2(node, edge, [], [], hfun);
##
##    figure;
##    patch('faces', tria(:,1:3), 'vertices', vert, 'facecolor', 'w', 'edgecolor', [.2,.2,.2]);
##    hold on;
##    plot(node(:,1), node(:,2), 'ro', 'markerfacecolor', 'r');
##    title(['Refined Mesh: ', title_str]);
##    axis equal;
end

