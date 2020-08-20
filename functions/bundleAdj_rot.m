function [Points, Camera] = bundleAdj_rot(Points, Camera, Obs, K)

MAXITR = 50; %%BA itr
TOL = 1e-6; %Step Tolerance

PARAMS{1}=size(Points,1);
PARAMS{2}=size(Camera,3);
PARAMS{3}=Obs;
PARAMS{4}=K;


X0=reshape(Points',1,3*PARAMS{1});



%{
for i=1:size(Camera,3)
   % Fill with i-th camera parameters
   X0 = [X0, ...
         rotm2eul(Camera(1:3,1:3,i)), ...
         Camera(1:3,4,i)'];

end

X = X0;
npoints = length(Obs{1}); %number of points
J   = zeros(2*PARAMS{2}*npoints, 6*PARAMS{2}+3*npoints);

mu = 10;

for itr = 1:MAXITR
    %Jacobian Matrix, obs shaping for r(p)
    for i=1:PARAMS{2}
        for j=1:size(Obs{i},1)
            point = X((1:3) + 3*(Obs{i}(3,j)-1) )';
            %camparam = [K(1,1), K(2,2), K(1,3), K(2,3), X(3*PARAMS{1}+(1:3)+6*(i-1)), X(3*PARAMS{1}+(4:6)+6*(i-1))];
            camparam = [X(3*PARAMS{1}+(1:3)+6*(i-1)), X(3*PARAMS{1}+(4:6)+6*(i-1))];
            
            %J_A = computeJacobianCameraParam(point',camparam);
            %J_B = computeJacobianPoint(point',camparam);
            J_A = JacobianCam(point',camparam, K);
            J_B = JacobianPoint(point',camparam,K);
            
            J(2*(j-1)+1+2*npoints*(i-1):2*(j-1)+2+2*npoints*(i-1), 1+6*(i-1):6+6*(i-1)) = J_A;
            J(2*(j-1)+1+2*npoints*(i-1):2*(j-1)+2+2*npoints*(i-1), 12+3*(j-1)+1:12+3*(j-1)+3) = J_B;
            
            %J(4*(j-1)+1+2*(i-1):4*(j-1)+2+2*(i-1), 1+6*(i-1):6+6*(i-1)) = J_A;
            %J(4*(j-1)+1+2*(i-1):4*(j-1)+2+2*(i-1), 12+3*(j-1)+1:12+3*(j-1)+3) = J_B;  
        end
    end
    
    if itr == 1
        lambda = 1;
    end
    [prevcost, rp] = computeCost(X, PARAMS, 1);
    

    
    H = J'*J;
    %H = H + eye(length(H)).*1e-10;
    H_mu = H + lambda * diag(diag(H));
    %H_mu = H + lambda * eye(length(H));
    J_rp = J'*rp;
    
    Xc_rp = J_rp(1:PARAMS{2}*6);
    Xp_rp = J_rp(PARAMS{2}*6+1:end);
    
    B = H_mu(1:PARAMS{2}*6,1:PARAMS{2}*6);
    E = H_mu(1:PARAMS{2}*6,PARAMS{2}*6+1:end);
    C = H_mu(PARAMS{2}*6+1:end, PARAMS{2}*6+1:end);
    S = B-E*inv(C)*E';
    if det(S) == 0
        S_inv = adjoint(S)/1e-6;
    else
        S_inv = inv(S);
    end
    
    deltaXc = S_inv*(Xc_rp - E*inv(C)*Xp_rp);
    deltaXp = inv(C)*(Xp_rp - E'*deltaXc);
    
    delta = [deltaXp' deltaXc'];
    
    X_new = X - delta;
    
    [curcost, rp] = computeCost(X_new, PARAMS, 1);
    
    if curcost < prevcost
        if (prevcost-curcost) <= TOL
            break;
        else
            X = X_new;
        end
        lambda = lambda/mu;
    else
        while curcost >= prevcost
            lambda  = lambda*mu;
            H_mu = H + lambda * diag(diag(H));
            %H_mu = H + lambda * eye(length(H));
            J_rp = J'*rp;
    
            Xc_rp = J_rp(1:PARAMS{2}*6);
            Xp_rp = J_rp(PARAMS{2}*6+1:end);
    
            B = H_mu(1:PARAMS{2}*6,1:PARAMS{2}*6);
            E = H_mu(1:PARAMS{2}*6,PARAMS{2}*6+1:end);
            C = H_mu(PARAMS{2}*6+1:end, PARAMS{2}*6+1:end);
            S = B-E*inv(C)*E';
    
            deltaXc = inv(S)*(Xc_rp - E*inv(C)*Xp_rp);
            deltaXp = inv(C)*(Xp_rp - E'*deltaXc);
    
            delta = [deltaXp' deltaXc'];
    
            X_new = X - delta;
            [curcost, rp] = computeCost(X_new,PARAMS, 1);

            disp(['prevcost: ',num2str(prevcost),'   curcost ',num2str(curcost),'   lambda ',num2str(lambda)]);
        end
        X = X_new;
    end
    disp(['Iteration: ',num2str(itr),'   Cost value is ',num2str(prevcost)]);
    
end
X = X_new;

Points = reshape(X(1:3*PARAMS{1}),[3 PARAMS{1}]);
Points = Points';
X(1:3*PARAMS{1}) = [];
for i=1:size(Camera,3)
   Camera(1:3,1:3,i) = eul2rotm(X(1:3)); 
   Camera(1:3,4,i)   = X(4:6)';
   X(1:6)            = [];
end

%}
end