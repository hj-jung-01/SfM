function [cost, rp] = computeCost(X, PARAMS, opt)
%opt 1 for rotm2eul
%opt 2 for rotm2quat

Obs = PARAMS{3};
K = PARAMS{4};
npoints = length(Obs{1});


if opt == 1
count = 1;
for i=1:PARAMS{2}
    rotat = eul2rotm( X(3*PARAMS{1}+(1:3)+6*(i-1)) );
    trans = X(3*PARAMS{1}+(4:6)+6*(i-1))';
    for j=1:size(Obs{i},1)
        point = X((1:3) + 3*(Obs{i}(j,3)-1) )';
        err(count:count+1) = Obs{i}(j,1:2) - proj(point, [rotat trans], K)';
        %rp(4*(j-1)+1+2*(i-1):4*(j-1)+2+2*(i-1)) = err(count:count+1);
        rp(2*(j-1)+1+2*npoints*(i-1):2*(j-1)+2+2*npoints*(i-1)) = err(count:count+1);
        count = count+2;
    end
end
cost = norm(err,2)^2;
rp = rp';

end

if opt == 2
    
count = 1;
for i=1:PARAMS{2}
    rotat = quat2rotm( X(3*PARAMS{1}+(1:4)+6*(i-1)) );
    trans = X(3*PARAMS{1}+(5:7)+6*(i-1))';
    
    for j=1:size(Obs{i},1)
        point = X((1:3) + 3*(Obs{i}(j,3)-1) )';
        err(count:count+1) = Obs{i}(j,1:2) - proj(point, [rotat trans], K)';
        %rp(4*(j-1)+1+2*(i-1):4*(j-1)+2+2*(i-1)) = err(count:count+1);
        rp(2*(j-1)+1+2*npoints*(i-1):2*(j-1)+2+2*npoints*(i-1)) = err(count:count+1);
        count = count+2;
    end
end
cost = norm(err,2)^2;
rp = rp';

end

end