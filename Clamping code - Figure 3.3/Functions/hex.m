% hex.m Outputs vertices of a hexagagonal grid, defined by number of rows,
% columns and side length n of the hexagons.

function [X_vertices, Y_vertices]=hex(rows,cols,n)

x_offset = n*sqrt(3);
y_offset = n+n/2;

% determine the center locations
X_centers=repmat((1:cols),rows,1);
X_centers = X_centers*x_offset;

Y_centers = repmat((1:rows)',1,cols);
Y_centers = Y_centers*y_offset;

% now shift odd rows over
odd_offset=n*sqrt(3)/2;
X_centers(rows-1:-2:1,:)=X_centers(rows-1:-2:1,:)+odd_offset;

X_vertices = zeros(rows,cols,6);
Y_vertices = zeros(rows,cols,6);

% topleft
X_vertices(:,:,1) = X_centers-n*sqrt(3)/2;
Y_vertices(:,:,1) = Y_centers+n/2;

% top
X_vertices(:,:,2) = X_centers;
Y_vertices(:,:,2) = Y_centers+n;

% topright
X_vertices(:,:,3) = X_centers+n*sqrt(3)/2;
Y_vertices(:,:,3) = Y_centers+n/2;

% botright
X_vertices(:,:,4) = X_centers+n*sqrt(3)/2;
Y_vertices(:,:,4) = Y_centers-n/2;

% bot
X_vertices(:,:,5) = X_centers;
Y_vertices(:,:,5) = Y_centers-n;

% botleft
X_vertices(:,:,6) = X_centers-n*sqrt(3)/2;
Y_vertices(:,:,6) = Y_centers-n/2;

% reshape vertices
X_vertices=permute(X_vertices,[3 1 2]);
X_vertices=reshape(X_vertices,6,[]);

Y_vertices=permute(Y_vertices,[3 1 2]);
Y_vertices=reshape(Y_vertices,6,[]);

end