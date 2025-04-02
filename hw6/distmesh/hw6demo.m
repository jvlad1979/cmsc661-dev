function [ p, t ] = distmesh_custom_demos( cases, do_plot )
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

% Distance function for the L-shape
fd = @(p) l_dpolygon(p, L_vertices);

% Mesh density function (uniform here)
%fh = @(p) ones(size(p,1),1);
fh = @(p) 0.8 + 0.4 * sqrt(p(:,1).^2 + p(:,2).^2);

% Generate the mesh
[p,t] = distmesh(fd, fh, 0.01, [0,0;1,1], L_vertices);
if( do_plot )
  clf
  patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
  title(s)
  axis tight
  axis equal
end
end

function [p,t] = test_case2( do_plot )
% Pentagon inside a pentagon
s = 'Custom Example 2: Pentagon inside a pentagon';
disp(s)

% Outer pentagon vertices
outer_pv = [0.5+cos((0:4)'*2*pi/5 + pi/2), 0.5+sin((0:4)'*2*pi/5 + pi/2)];

% Inner pentagon vertices (scaled down)
inner_pv = 0.5 * ([0.5+cos((0:4)'*2*pi/5), 0.5+sin((0:4)'*2*pi/5)] - 0.5) + 0.5;

% Custom distance function for nested pentagons
fd = @(p) max( l_dpolygon(p, outer_pv), -l_dpolygon(p, inner_pv) );

fh = @(p) ones(size(p,1),1);
[p,t] = distmesh( fd, fh, 0.05, [0,0;1,1], [outer_pv; inner_pv] );

if( do_plot )
  clf
  patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
  title(s)
  axis tight
  axis equal
end
end

function [p,t] = test_case3( do_plot )
% Circle with interior circle removed
s = 'Custom Example 2: Circle with interior circle';
disp(s)

% Define the distance function for a circle with an interior circle removed
% Outer circle radius 1, centered at (0.5, 0.5)
% Inner circle radius 0.3, also centered at (0.5, 0.5)
fd = @(p) max(max( max(sqrt((p(:,1)-0.5).^2 + (p(:,2)-0.5).^2) - 1,p(:,2)-.5), ... % Outer circle
               -(sqrt((p(:,1)-0.25).^2 + (p(:,2)-0.2).^2) - 0.2)
               ),-(sqrt((p(:,1)-0.75).^2 + (p(:,2)-0.2).^2) - 0.2)); % Inner circle

% Uniform mesh size function
fh = @(p) ones(size(p,1),1);

% Generate mesh using distmesh
[p,t] = distmesh( fd, fh, 0.01, [0,0;1,1], [] );

% Optional plotting
if( do_plot )
  clf
  patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )
  title(s)
  axis tight
  axis equal
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
