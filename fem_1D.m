clear 
close all
clc


% INPUT DATA
l = 1;      % length of segment
h = 0.01;   % step
F = 1;      % force applied homogeneously
m = 10;     % mass of the object
% BOUNDARY CONDITIONS
alpha = 0;  % x = 0
beta = 0;   % x = l
% 1 : Dirichlet, 2 : Neumann 
type0 = 1;
type1 = 2;

% Main code
x = 0:h:l;
nodi = x(1):h:x(end);
N = length(nodi);
fun = -F*ones(1,N);
A = (m/h^2)*((diag(2*ones(1,N)) + diag(-1*ones(1,N-1),-1) ...
    + diag(-1*ones(1,N-1),1)));
% BC
if type0 == 1
    A(1,:) = zeros(1,N);
    A(1,1) = 1;
elseif type0 == 2
    A(1,1:3) = (1/(2*h))*[-3 4 -1];
end

if type1 == 1
    A(end,:) = zeros(1,N);
    A(end,end) = 1;
elseif type1 == 2
    A(end,end-2:end) = (1/(2*h))*[1 -4 3];
end

% Solution Vector
b = zeros(N,1);
b(2:N-1) = (fun(2:end-1));
b(1) = alpha;
b(end) = beta;

% Solving and plotting
u = A\b;
plot(x,u); axis equal;