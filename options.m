%options = optimoptions('ga','MutationFcn',{@myfun, rate})
function mutationChildren = myfun(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
   % SET G = 100 IN THIS CASE
   G=100;
  % mutation  
  r = rand(1,1);
  b = rand(1,1);
  mutationChildren=parents.*r;
  for k=1:nvars
    mutationChild1=parents(:,k)+(upper(parents(:,k))-parents(:,k))*r*(1-Generation/G)^b;
    mutationChild1=parents(:,k)+(upper(parents(:,k))-parents(:,k))*r*(1-Generation/G)^b;
    
    rndom = rand(1,1);
    if rndom>.5
      rndom=1;
    else 
      rndom=0;
    end
    
    mutationChildren(:,k) = rndom*mutationChild1+(1-rndom)*mutationChild2
  end
  
  % state structure
  rate = exp(-Generation/(4*G))-1;
  
end
%-------------------------------------------------------------------------------------------------
function y = lower(k)  
% which -all worksopacefunc
    if k<(N+1)
      y=-10;
    elseif k==(4*N+1)
      y=6.22e-15;
    else 
      y=-20;
    end
% which -all worksopacefunc
end

function y = upper(k)   
% which -all worksopacefunc
    if k<(N+1)
    y=40;
    elseif k==(4*N+1)
    y=10;
    else 
    y=20;
    end
% which -all worksopacefunc
end
%---------------------------------------------------------------------------------------
%[x,fval] = ga(@FitnessFcn,45,[],[],[],[],[],[],[],options);
