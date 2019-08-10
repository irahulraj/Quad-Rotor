function [child] = mutation(child)

    r = randn(1,1);
    
    if r<0 && r>1
    r = 0.599;
    end
 
    Gene_no = length(child.Gene);
    for k = 1: Gene_no
        R = rand();
        xk_dash1 = child(k) + (upper(x(k))- child(k))*r*(1-gen/G)^b; 
        xk_dash2 = child(k) + ((x(K)- lower(child(k))))*r*(1-gen/G)^b; 
        if R < Pm
            rn = randn(1,1);
                if rn<0.5
                rn = 0;
                else
                rn=1;
                end
            child.Gene(k) = xk_dash1*rn+(1-rn)*xk_dash2;
        end
    end


end