function [x1_dash, x2_dash] = crossover(x1,x2)
    lambda = randn(1,1);
    
    if lambda<0 && lambda>1
    lambda = 0.499;
    end

    x1_dash = lambda*x1+(1-lambda)*x2;
    x2_dash = lambda*x2+(1-lambda)*x1;
% disp(x1_dash)
% disp(x2_dash);
end
