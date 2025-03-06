% ras: compute coefficients and region of absolute stability for 
% classical linear multistep methods
% a_1 u_{n+1} + a_2 u_n + ... + a_{k+1} u_{n-k+1}
%  = h ( b_1 f_{n+1} + b_2 f_n + ... + b_{k+1} f_{n-k+1} )
%
% k-step explicit Adams of order p = k
% where a_1 =1, a_2 = -1, a_3 = ... = a_{k+1} = 0, b_1 = 0.
%
% k-step implicit Adams of order p = k + 1
% where a_1 =1, a_2 = -1, a_3 = ... = a_{k+1} = 0.
%
% k-step implicit BDF of order p = k
% where b_1 = 1, b_2 = ... = b_{k+1} = 0.
%
% Impose p + 1 order conditions rho(e^h) - h sigma(e^h) = O(h^{p+1})
% on the degree-k polynomials
% rho(z) = a_1 z^k + ... + a_{k+1} 
% sigma(z) = b_1 z^k + ... + b_{k+1} 
% plus the above k + 2 or k + 1 constraints on the a or b coeffs: then
% solve r = p + 1 + k + 2 (explicit or BDF) or p + 1 + k + 1 (implicit) 
% equations in c = 2k + 2 unknowns by least squares i.e. SVD.

I = complex(0,1);  % Set up an imaginary unit (Octave's I = Matlab's i)

% Loop through method m, step number k and order p.
% Method is 1 for explicit Adams, 2 for implicit Adams, 3 for BDF, ...
ead = 1;
iad = 2;
bdf = 3;
% File names for plot files e-2.ps and so forth.
cap = 'eib';

run( ead , 1 ) = ead;		% Method
run( ead , 2 ) = 9;					% Number of orders
run( ead , 3 ) = 0;					% p minus k

run( iad , 1 ) = iad;		% Method
run( iad , 2 ) = 9;					% Number of orders
run( iad , 3 ) = 1;					% p minus k

run( bdf , 1 ) = bdf;					% Method
run( bdf , 2 ) = 9;					% Number of orders
run( bdf , 3 ) = 0;					% p minus k

% Run through methods.
for test = 1 : 3;
% Run orders k = 2 up to maximum specified.
for k = 2 : run( test, 2 )

method = run( test, 1 )
p = k + run( test, 3 )

% Number of coeffs in a and in b.
n = k + 1;

% Size of matrix m of order conditions: r rows by c columns.
if method == ead
 r = p + 1 + k + 2;
 c = 2 * n;
end

if method == iad | method == bdf 
 r = p + 1 + k + 1;
 c = 2 * n;
end

m = zeros( r, c );

% First half-row for the a coeffs.
for j = 1 : n
 m( 1, j ) = 1;
end

% Remaining half-rows for the a coeffs.
for q = 2 : p + 1

 for j = 1 : n
  m( q , j ) = ( -j )^( q - 1 ) / gamma( q );
 end

end

% Second halves of first p + 1 rows are minus first halves shifted.
for jp = n + 1 : 2 * n
 j = jp - n;
 m( 1, jp ) = 0;
 for q = 2 : p + 1
  m( q, jp ) = - m( q - 1, j );
 end
end

% Next row to add to r by c matrix
next = p + 2;

if method == ead
% Explicit Adams: set b1 = 0, a1 = 1, a2 = -1, other a's equal to 0.

 m( next : r - 1, : ) = eye( r - next, c );

% Set b1 = 0 and print final m.
 m( r, n + 1 ) = 1

% Right-hand side of equations mostly 0.
rhs = zeros(r, 1);
rhs( next, 1 ) = 1;
rhs( next + 1, 1 ) = -1;
rhs( r, 1 ) = 0

end

if method == iad
% Implicit Adams: set a1 = 1, a2 = -1, other a's equal to 0.

 m( next : r, : ) = eye( r - next + 1, c )

% Right-hand side of equations mostly 0.
rhs = zeros(r, 1);
rhs( next, 1 ) = 1;
rhs( next + 1, 1 ) = -1;

end

% BDF and its variants: set b1 = 1, other b's equal to 0.
if method == bdf 

 m( next : r, n + 1 : c ) = eye( r - next + 1, c - n )

% Right-hand side of equations mostly 0.
rhs = zeros(r, 1);
rhs( n + 1, 1 ) = 1

end

% Compute least-squares solution x.
x = m \ rhs;

% Extract left and right coeffs of method.
a = x( 1 : n )
b = x( n + 1 : 2 * n )

% Store boundary locus in f for plotting

z = linspace(0, 2 * pi, 1000 );
z = exp( z * I );

av = polyval( a, z );
bv = polyval( b, z );
f = av ./ bv;
xmax = max( real( f ) );
xmin = min( real( f ) );
ymax = max( imag( f ) );
ymin = min( imag( f ) );

% Store grid points hL where all roots are inside unit circle in g 
g( 1 : 2500 ) = 0;
n = 1;
for ix = 0 : 50
 x = xmin + ix * ( xmax - xmin ) / 50.0;
 for iy = 0 : 50
  y = ymin + iy * ( ymax - ymin ) / 50.0;
  hL = x + I * y;
  stab = a - hL * b;
  root = roots( stab );
  absroot = abs( root );
  maxroot = max( absroot );
  if maxroot <= 1.0
   g( n ) = hL;
   n = n+1;
  end
 end
end

% Produce Fig output file: convert to PostScript by fig2dev -Leps

% gset term fig
% s = sprintf( "gset output '%s-%d';", cap( :, method ), k );
% eval( s )

% Plot boundary-including locus.
plot( f )
% Add points where roots are inside unit circle.
hold on
plot( g, '+' )
hold off
print('-dps', sprintf('%s-%d', cap(:, method), k) );

end
end