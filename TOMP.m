Ix = 0.0142;
Iy = 0.0142;
Iz = 0.0071;
m = 0.56;
l = 0.21;
C = 1.3;
del_t = 0.025;
N = 11;

% Declaring variables
phi_dot_dot = zeros(1,N);
phi_dot = zeros(1,N);
phi = zeros(1,N);
theta_dot_dot = zeros(1,N);
theta_dot = zeros(1,N);
theta = zeros(1,N);
psi_dot_dot = zeros(1,N);
psi_dot = zeros(1,N);
psi = zeros(1,N);

x_dot_dot = zeros(1,N); 
x_dot = zeros(1,N);
x = zeros(1,N);
y_dot_dot = zeros(1,N);
y_dot = zeros(1,N);
y = zeros(1,N);
z_dot_dot = zeros(1,N);
z_dot = zeros(1,N);
z = zeros(1,N);

% Initialization the variables
    % initial point
    phi_dot_dot(1) = 0;
    phi_dot(1) = 0;
    phi(1) = 0;
    theta_dot_dot(1) = 0;
    theta_dot(1) = 0;
    theta(1) = 0;
    psi_dot_dot(1) = 0;
    psi_dot(1) = 0;
    psi(1) = 0;
    
    x_dot_dot(1) = 0;
    x_dot(1) = 0;
    x(1) = 0.3;
    y_dot_dot(1) = 0;
    y_dot(1) = 0;
    y(1) = 0.4;
    z_dot_dot(1) = 0;
    z_dot(1) = 0;
    z(1) = 0.5;
    
    %final point
   % phi_dot(N) = 0;
   % phi(N) = 0;
    % theta_dot(N) = 0;
    %theta(N) = 0;
    %psi_dot(N) = 0;
    %psi(N) = 0;

    %x_dot(N) = 0;
%    x(N) = 0;
 %   y_dot(N) = 0;
  %  y(N) = 0;
   % z_dot(N) = 0;
    %z(N) = 0;

%-------------------------------------------------------------------------------
for i=2:N
    phi_dot_dot(i) = u3(i)*l/Ix; 
    phi_dot(i) = phi_dot(i-1) + phi_dot_dot(i-1)*del_t;
    phi(i) = phi(i-1) + phi_dot(i-1)*del_t + (phi_dot_dot(i-1))*del_t^2/2;

    theta_dot_dot(i) = u2(i)*l/Iy; 
    theta_dot(i) = theta_dot(i-1) + theta_dot_dot(i-1)*del_t;
    theta(i) = theta(i-1) + theta_dot(i-1)*del_t + (theta_dot_dot(i-1))*del_t^2/2;

    psi_dot_dot(i) = u4(i)*C/Iz; 
    psi_dot(i) = psi_dot(i-1) + psi_dot_dot(i-1)*del_t;
    psi(i) = psi(i-1) + psi_dot(i-1)*del_t + (psi_dot_dot(i-1))*del_t^2/2;


    x_dot_dot(i) = (sin(psi(i))*sin(phi(i))+cos(psi(i))*sin(theta(i))*cos(phi(i)))*u1(i)/m; 
    x_dot(i) = x_dot(i-1) + x_dot_dot(i-1)*del_t;
    x(i) = x(i-1) + x_dot(i-1)*del_t + (x_dot_dot(i-1))*del_t^2/2;

    y_dot_dot(i) = (sin(psi(i))*sin(phi(i))+cos(psi(i))*sin(theta(i))*cos(phi(i)))*u1(i)/m; 
    y_dot(i) = y_dot(i-1) + y_dot_dot(i-1)*del_t;
    y(i) = y(i-1) + y_dot(i-1)*del_t + (y_dot_dot(i-1))*del_t^2/2;

    z_dot_dot(i) = (sin(psi(i))*sin(phi(i))+cos(psi(i))*sin(theta(i))*cos(phi(i)))*u1(i)/m; 
    z_dot(i) = z_dot(i-1) + z_dot_dot(i-1)*del_t;
    z(i) = z(i-1) + z_dot(i-1)*del_t + (z_dot_dot(i-1))*del_t^2/2;
end   


x_f=0;
y_f=0;
z_f=0;
phi_f=0;
theta_f=0;
psi_f=0;