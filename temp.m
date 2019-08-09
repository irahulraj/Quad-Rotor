% ------------------------------------------------------------------------------
phi_dot_dot(i) = u3(i)*l/Ix; 
phi_dot(i) = phi_dot(1) + (l/Ix)*(sum(u3(1:i-1));
j=1:i-1;
phi(i) = phi(1) + (i-1)*phi_dot(1)*del_t+sum(sum(u3(1:j-1)*l/Ix))*del_t^2+(sum(u3(1:j-1)*l/Ix))/2;

theta_dot_dot(i) = u2(i)*l/Iy; 
theta_dot(i) = theta_dot(1) + (l/Iy)*(sum(u2(1:i-1));
j=1:i-1;
theta(i) = theta(i-1) + (i-1)*theta_dot(1)*del_t+sum(sum(u2(1:j-1)*l/Iy))*del_t^2+(sum(u2(1:j-1)*l/Iy))/2;

psi_dot_dot(i) = u4(i)*C/Iz; 
psi_dot(i) = psi_dot(1) + (C/Iz)*(sum(u3(1:i-1));
j=1:i-1;
psi(i) = psi(i-1) + (i-1)*psi_dot(1)*del_t+sum(sum(u3(1:j-1)*C/Iz))*del_t^2+(sum(u3(1:j-1)*C/Iz))/2;


x_dot_dot(i) = (sin(psi(i))*sin(phi(i))+cos(psi(i))*sin(theta(i))*cos(phi(i)))/m; 
x_dot(i) = phi_dot(1) + (l/Ix)*(sum(u3(1:i-1));
j=1:i-1;
x(i) = phi(1) + (i-1)*x_dot(1)*del_t+sum(sum(u3(1:j-1)*l/Ix))*del_t^2+(sum(u3(1:j-1)*l/Ix))/2;

% -------------------------------------------------------------------------------------------------------
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
