function [ p, t ] = hw6( cases, do_plot )
% DISTMESH_CUSTOM_DEMOS Run custom DistMesh example demos.
%
%   Examples:
%   1. L-shaped region
%   2. Pentagon inside a pentagon
%   3. Semicircle with two holes

close all
clc

if( ~nargin || isempty(cases) )
  cases = 1:3;
end
if( nargin<2 )
  do_plot = true;
end

for i=cases
  eval( ['[p,t] = test_case',num2str(i),'(',num2str(do_plot),');'] );
  if( do_plot && i~=cases(end) )
    pause
  end
  fprintf('\n\n')
end

if( ~nargout )
  clear p t
end
end

function [p,t] = test_case1( do_plot )
% L-shaped region
s = 'Custom Example 1: L-shaped region';
disp(s)

% Define L-shape using vertices
L_vertices = [0 0; 1 0; 1 0.4; 0.4 0.4; 0.4 1; 0 1]; % L shape
L_edges = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1]; % Connectivity list

% Distance function for the L-shape
fd = @(p) l_dpolygon(p, L_vertices);

% Mesh density function (uniform here)
%fh = @(p) ones(size(p,1),1);
fh = @(p) 0.8 + 0.4 * sqrt(p(:,1).^2 + p(:,2).^2);

% Generate the mesh
[p,t] = distmesh(fd, fh, 0.01, [0,0;1,1], L_vertices);

node = L_vertices
edge = L_edges

##edges = t
[vert, etri, tria, tnum] = refine2(node, edge);



if do_plot
    figure(1); % First figure
    clf; % Clear only this figure
    patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9]);
    title(s);
    axis tight;
    axis equal;
end

[p, t] = mesh2d(node, edge);
if do_plot
    figure(2); % First figure
    clf; % Clear only this figure
    patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9]);
    title(s);
    axis tight;
    axis equal;
end



end

function [p,t] = test_case2( do_plot )
% Pentagon inside a pentagon
s = 'Custom Example 2: Pentagon inside a pentagon';
disp(s)

% Outer pentagon vertices
outer_pv = [0.5+cos((0:4)'*2*pi/5 + pi/2), 0.5+sin((0:4)'*2*pi/5 + pi/2)];
num_vertices = size(outer_pv, 1);
% Inner pentagon vertices (scaled down)
inner_pv = 0.5 * ([0.5+cos((0:4)'*2*pi/5), 0.5+sin((0:4)'*2*pi/5)] - 0.5) + 0.5;
num_vertices = size(outer_pv, 1);
outer_edges = [(1:num_vertices)' [2:num_vertices 1]'];  % Outer pentagon
inner_edges = [(1:num_vertices)' [2:num_vertices 1]'];  % Inner pentagon

node = [outer_pv; inner_pv]
edge = [outer_edges;inner_edges+num_vertices]

% Custom distance function for nested pentagons
fd = @(p) max( l_dpolygon(p, outer_pv), -l_dpolygon(p, inner_pv) );

fh = @(p) ones(size(p,1),1);
[p,t] = distmesh( fd, fh, 0.05, [0,0;1,1], [outer_pv; inner_pv] );



##[vert, etri, tria, tnum] = refine2(node, edge);

if do_plot
    figure(1); % First figure
    clf; % Clear only this figure
    patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9]);
    title(s);
    axis tight;
    axis equal;
end

[p, t] = mesh2d(node, edge);
if do_plot
    figure(2); % First figure
    clf; % Clear only this figure
    patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9]);
    title(s);
    axis tight;
    axis equal;
end
end
function [p, t] = test_case3(do_plot)
    % Circle with interior cut-outs
    s = 'Custom Example 2: Circle with interior cut-outs';
    disp(s);

    % Outer boundary (semicircle + straight line at the top)
    num_outer = 64;   % Points on the outer semicircle
    num_inner = 16;   % Points on each inner cut-out

    % Outer semicircle
    theta_outer = linspace(pi, 2*pi, num_outer)';
    outer_nodes = [0.5 + cos(theta_outer), 0.5 + sin(theta_outer)];

    % Add straight-line segment at the top to close the semicircle
    top_line = [outer_nodes(1,:); outer_nodes(end,:)];

    % Inner circles (centered at (0.25,0.2) and (0.75,0.2), radius = 0.2)
    theta_inner = linspace(0, 2*pi, num_inner)';
    inner1_nodes = [0.25 + 0.2 * cos(theta_inner), 0.2 + 0.2 * sin(theta_inner)];
    inner2_nodes = [0.75 + 0.2 * cos(theta_inner), 0.2 + 0.2 * sin(theta_inner)];

    % Combine all nodes with separate indexing
    node = [
        outer_nodes;     % 1 to num_outer
        top_line;        % num_outer+1 to num_outer+2
        inner1_nodes;    % num_outer+3 to num_outer+num_inner+2
        inner2_nodes     % num_outer+num_inner+3 to num_outer+2*num_inner+2
    ];

    % Create edges with corrected indexing
    % Outer semicircle edges
    outer_edges = [
        (1:num_outer-1)', (2:num_outer)';       % Semicircle curve
        num_outer, num_outer+1;                 % To top line
        num_outer+1, 1                          % Close the outer shape
    ];

    % Inner circle 1 edges
    inner1_edges = [
        (num_outer+3:num_outer+num_inner+2-1)', (num_outer+4:num_outer+num_inner+2)';  % Inner circle 1
        num_outer+num_inner+2, num_outer+3      % Close inner circle 1
    ];

    % Inner circle 2 edges
    inner2_edges = [
        (num_outer+num_inner+3:num_outer+2*num_inner+2-1)', (num_outer+num_inner+4:num_outer+2*num_inner+2)';  % Inner circle 2
        num_outer+2*num_inner+2, num_outer+num_inner+3  % Close inner circle 2
    ];

    % Combine all edges
    edge = [outer_edges; inner1_edges; inner2_edges];

    % Signed distance function using l_dpolygon for each separate polygon
    fd_outer = @(p) l_dpolygon(p, [outer_nodes; flipud(top_line)]); % Ensure correct ordering
    fd_inner1 = @(p) -l_dpolygon(p, inner1_nodes); % Negative for holes
    fd_inner2 = @(p) -l_dpolygon(p, inner2_nodes); % Negative for holes

    % Final distance function (union of boundaries)
    fd = @(p) max(fd_outer(p), max(fd_inner1(p), fd_inner2(p)));

    % Uniform mesh size function
    fh = @(p) ones(size(p,1),1) * 0.1;

    % Generate mesh using distmesh
    [p, t] = distmesh(fd, fh, 0.05, [0, 0; 1, 1], []);

    part = {
        1:size(outer_edges,1),                     % Outer semicircle region
        size(outer_edges,1)+1:size(outer_edges,1)+size(inner1_edges,1),  % First inner circle
        size(outer_edges,1)+size(inner1_edges,1)+1:size(edge,1)  % Second inner circle
    };
    figure;
    hold on;
    plot(node(:,1), node(:,2), 'bo'); % Plot nodes
    for k = 1:size(edge,1)
      plot(node(edge(k,:),1), node(edge(k,:),2), 'r-');
    end
    hold off;

    % Optional: Refinement options
    opts.kind = 'DELFRONT';   % Frontal-Delaunay refinement
    opts.rho2 = 1.025;        % Radius-edge ratio threshold
    opts.siz2 = 1.300;        % Relative length threshold for triangles

##    [vert, etri, tria, tnum] = refine2(inner1_nodes, inner1_edges);

    if do_plot
        figure(1); % First figure
        clf; % Clear only this figure
        patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9]);
        title(s);
        axis tight;
        axis equal;
    end
    [p, t] = mesh2d(node, edge);
    if do_plot
        figure(2); % First figure
        clf; % Clear only this figure
        patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9]);
        title(s);
        axis tight;
        axis equal;
    end
end


% Import the existing l_dpolygon function from the original script
function [ dist ] = l_dpolygon( p, v )
    n_p = size(p,1);
    n_s = size(v,1)-1;

    dist = zeros(n_s,n_p);
    for i_s=1:n_s
      v_i = v([i_s,i_s+1],:);
      n_p = size(p,1);
      w  = v_i(2,:)-v_i(1,:);
      ix1 = ones(n_p,1);
      vp = v_i(ix1,:)-p;
      w1 = w(ix1,:);
      s = dot(w1,vp,2);
      u = -s/(w*w.');
      u(u<0) = 0;
      u(u>1) = 1;
      h = w1.*[u,u]+vp;
      dist(i_s,:) = sqrt(dot(h,h,2));
    end
    dist = (-1).^(inpolygon(p(:,1),p(:,2),v(:,1),v(:,2))).*min(dist).';
end
