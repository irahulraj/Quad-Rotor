% Bounds
N = 11;
GA_LB1 = -10*ones(1,N);
GA_LB2 = -20*ones(1,N);
GA_LB3 = -20*ones(1,N);
GA_LB4 = -20*ones(1,N);

GA_UB1 = 40*ones(1,N);
GA_UB2 = 20*ones(1,N);
GA_UB3 = 20*ones(1,N);
GA_UB4 = 20*ones(1,N);

GA_LB = [GA_LB1, GA_LB2, GA_LB3, GA_LB4, 6.22e-15];
GA_UB = [GA_UB1, GA_UB2, GA_UB3, GA_UB4, 10];
