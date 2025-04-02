rand('state',1); % Always the same results
set(gcf, 'renderer', 'opengl');
%fprintf('Rectangle with circular hole, refined at circle boundary\n');
%fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
%fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
%[p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%fprintf('(press any key)\n\n'); pause

pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
    1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
[p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);
