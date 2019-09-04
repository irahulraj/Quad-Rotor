% GA algorithms

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


% fitness function 
%function Y = fit_func()
%e_sqr = ((x_f-x(N))^2+(y_f-y(N))^2+(z_f-z(N))^2+(phi_f-phi(N))^2+(theta_f-theta(N))^2+(psi_f-psi(N))^2);
%e_dot_sqr = ((x_dot(N))^2+(y_dot(N))^2+(z_dot(N))^2+(phi_dot(N))^2+(theta_dot(N))^2+(psi_dot(N))^2);
%Y = -1/(1+e_sqr+e_dot_sqr);
%end
function fitness = fit_func(chromosome)
  e_sqr = ((x_f-x(N))^2+(y_f-y(N))^2+(z_f-z(N))^2+(phi_f-phi(N))^2+(theta_f-theta(N))^2+(psi_f-psi(N))^2);
  e_dot_sqr = ((x_dot(N))^2+(y_dot(N))^2+(z_dot(N))^2+(phi_dot(N))^2+(theta_dot(N))^2+(psi_dot(N))^2);
  fitness = 1/(1+e_sqr+e_dot_sqr);
endfunction
pop_size = 50;


% Produce Initial generation in random way
% -------HERE----

% OffSpring generation

% CrossOver
function [child1 , child2] = crossover(parent1 , parent2, Pc)
%  Using whole arithmetic crossover  ref - http://webpages.iust.ac.ir/yaghini/Courses/AOR_872/Genetic%20Algorithms_03.pdf
%  lambda = rand();  
  child1 = lambda*parent1 + (1-lambda)*parent2;
  child2 = lambda*parent2 + (1-lambda)*parent1;  
  
  R1 = rand();
  if R1 <= Pc
    child1 = child1;
  else
    child1 = parent1;
  endif
  
  R2 = rand();
  if R2 <= Pc
    child2 = child2;
  else
    child2 = parent2;
  endif
  % 
  
endfunction
%

% Mutation
function child = mutation(parent, Pm)
%  r = rand();
  Gene_no = 4*N+1;
  for k = 1:Gene_no
%    R = rand();
    if(R>Pm)
      child1(k) = parent(k) + (upper_lim(k) - parent(k))*r*(1-gen/G)^b;
      child2(k) = parent(k) - (parent(k) - lower_lim(k))*r*(1-gen/G)^b;
    else
      continue;
    endif
  endfor
    %
  %  selection = rand();
  if(selection>0.5)
    child = child1;
  else 
    child = child2;
  endif
  %  
endfunction 
%


% Selectin methods  - Rank Mechanism
%https://stackoverflow.com/questions/34961489/rank-selection-in-ga
NewFitness=sort(Fitness);
        NewPop=round(rand(PopLength,IndLength));

        for i=1:PopLength
            for j=1:PopLength
                if(NewFitness(i)==Fitness(j))
                    NewPop(i,1:IndLength)=CurrentPop(j,1:IndLength);
                    break;
                end
            end
        end
        CurrentPop=NewPop;

        ProbSelection=zeros(PopLength,1);
        CumProb=zeros(PopLength,1);

        for i=1:PopLength
            ProbSelection(i)=i/PopLength;
            if i==1
                CumProb(i)=ProbSelection(i);
            else
                CumProb(i)=CumProb(i-1)+ProbSelection(i);
            end
        end

        SelectInd=rand(PopLength,1);

        for i=1:PopLength
            flag=0;
            for j=1:PopLength
                if(CumProb(j)<SelectInd(i) && CumProb(j+1)>=SelectInd(i))
%                    SelectedPop(i,1:IndLength)=CurrentPop(j+1,1:IndLength);
                    selected_parent1 = CurrentPop(j+1,1:IndLength);
                    flag=1;
                    break;
                end
            end
            if(flag==0)
                SelectedPop(i,1:IndLength)=CurrentPop(1,1:IndLength);
            end
        end
%










% main function

% Random Population Generation
rows = pop_size;
columns = 4*N+1;
Numbers = rand(rows,columns)

Numbers(:,1:N) = Numbers(:N)*50-10;
Numbers(:,N+1:2*N) = Numbers(:N)*50-10;
Numbers(:,2*N+1:3*N) = Numbers(:N)*50-10;
Numbers(:,3*N+1:4*N) = Numbers(:N)*50-10;

fitness_val = zeros(pop_size,1);

for i = 1:pop_size
  fitness_val(k) = fitness(population(k));
  endfor
%
% maximal generation
  G = 100;

for gen = 1:G
  crossover_rate = exp(-gen/G);
  mutation_rate = exp(-gen/(4*G))-1;
  
  new_population = zeros(pop_size, Gene_no);
  % Crossover 
  for pop = 1:pop_size
    new_population(pop:pop+1,:) = crossover(population(pop), population(pop+1), crossover_rate);
  endfor
  %
  for pop = 1:pop_size
    new_population(pop) = mutation(population(pop), mutation_rate);
  endfor
  
  
  
    
  
  
  endfor
%
  























