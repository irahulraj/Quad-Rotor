function Y = fit_func(XU)

% XU = [u11,....u1N,u21,......,u31,.....,u41,....., del_t] = 4N+1 vars
% constraints =
% GA_LB = [-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,-20.0,6.2e-11]
% GA_UB = [40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,9999999999999999.0]
    
%      give N through command window
    N = 11;
    bounds;
    u1 = XU(1:N);
    u2 = XU(N+1:2*N);
    u3 = XU(2*N+1:3*N);
    u4 = XU(3*N+1:4*N);    
    del_t = XU(4*N+1);
    
    Ix = 0.0142;
    Iy = 0.0142;
    Iz = 0.0071;
    m = 0.56;
    l = 0.21;
    C = 1.3;
%     del_t = 0.025;
    

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
    
    
e_sqr = ((x_f-x(N))^2+(y_f-y(N))^2+(z_f-z(N))^2+(phi_f-phi(N))^2+(theta_f-theta(N))^2+(psi_f-psi(N))^2);
e_dot_sqr = ((x_dot(N))^2+(y_dot(N))^2+(z_dot(N))^2+(phi_dot(N))^2+(theta_dot(N))^2+(psi_dot(N))^2);
Y = -1/(1+e_sqr+e_dot_sqr);
end