% Mesh for bump in channel

clear all
close all
clc

tic
l_x = 1.5; % length in X
l_y = 1.5; % length in Y
n_x = 36;  % number of points in X
n_y = 36;  % number of points in Y

x = linspace(0, l_x, n_x+1);
dx = x(2) - x(1); % dx

y = zeros(1, n_x+1);

b_1 = int16(0.3/dx);
b_2 = int16(1.2/dx);

for i = b_1:b_2
    y(i+1) = 0.05*(sin(pi*x(i+1)/0.9 - (pi/3.)))^4;%comment this line for a flat mesh
end

d_k = l_y/2;  % Distance from the viscous wall where the mesh will be denser
p_k = 50;     % Percentage concentration for the zone with higher mesh density
n_yp = n_y*(p_k/100); % Number of nodes in the high-density region
n_yf = n_y - n_yp;    % Number of nodes in the free region

dy_p = d_k / n_yp;     % y-step in the high-density part of the mesh
dy_f = (l_y - d_k) / n_yf; % y-step in the remaining part of the mesh

% Matrix that will store the X coordinates of the mesh
X = zeros(n_y+1, n_x+1);
% Storing X coordinates, constant step
for i = 1:n_x+1
    X(:, i) = x(i);
end

% Initial matrix that will store the Y coordinates of the mesh
Y_t = zeros(n_y+1, n_x+1);
for i = 1:n_x+1
    Y_t(1, i) = y(i);
end

% Storing Y coordinates for the high-density nodes
for i = 2:n_yp+1
    for j = 1:n_x+1
        Y_t(i, j) = Y_t(i-1, j) + dy_p;
    end
end

% Storing Y coordinates for the rest of the mesh, dynamic dy_f
dy_f_d = 0;          % dynamic y-step
max_y = max(y);      % maximum value of Y
dy_din = zeros(1, n_x+1); % vector for y-steps in the free part
for i = n_yp+2:n_y+1
    for j = 1:n_x+1
        % dy_f_d = (l_y + max_y - y(j)) / n_yf;
        dy_f_d = (l_y - (y(j) + d_k)) / n_yf;
        Y_t(i, j) = Y_t(i-1, j) + dy_f_d;
        dy_din(j) = dy_f_d / dx;
    end
end

Y = zeros(size(Y_t));

k = 0;
for i = 1:n_y+1
    for j = 1:n_x+1
        Y(i, j) = Y_t(n_y+1-k, j);
    end
    k = k + 1;
end

clear Y_t;

z = zeros(n_y+1, n_x+1);

figure(1)
plot(x, y)

% Flatten the matrices into vectors
x_coords = X(:);
y_coords = Y(:);
%{
% Create a scatter plot
figure(2)
scatter(x_coords, y_coords, 10, 'filled');  % The '10' specifies the marker size
% Label the axes
xlabel('X Coordinates');
ylabel('Y Coordinates');
% Title for the plot
title('Scatter Plot of Unstructured Mesh Points');
% Set axis equal for better representation of the points
axis equal;
% Optionally, add grid and adjust other plot settings
grid on;
%}
figure(3)
mesh(X, Y, z)
title('Mesh')
axis equal

%{
% Derivative of Y
y_p = 0;
for i = 1:n_x
    y_p(i) = (y(i+1) - y(i)) / dx;
end
figure(2)
plot(x(1:n_x), y_p)
%}

% dy/dx relation in the free part
figure(4)
plot(1:n_x+1, dy_din)

% Volume of the cells
vol_f = zeros(1, n_x);
vol_b = zeros(1, n_x);
for i = 1:n_x
    vol_f(i) = -0.5 * ((Y(n_y+1, i) - Y(n_y, i)) + (Y(n_y+1, i+1) - Y(n_y, i+1))) * dx;
    vol_b(i) = -0.5 * ((Y(2, i) - Y(1, i)) + (Y(2, i+1) - Y(1, i+1))) * dx;
end

figure(5)
plot(1:n_x, vol_f)
hold on 
plot(1:n_x, vol_b)
hold off
title("Volumes of the cells")

% C: distance in X to the cell's centroid
c_cent = zeros(1, n_x);
for i = 1:n_x
    a_i = Y(1, i) - Y(2, i);
    b_i = Y(1, i+1) - Y(2, i+1);
    c_cent(i) = a_i / (b_i + a_i);
end

figure(6)
plot(1:n_x, c_cent)
title("Relation C/dx")

save('mesh_bumpchannel2.mat', 'X', 'Y');

toc
