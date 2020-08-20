function [Points, Camera] = bundleAdj_quat(Points, Camera, Obs, K)


PARAMS{1}=size(Points,2);
PARAMS{2}=size(Camera,1);
PARAMS{3}=Obs;
PARAMS{4}=K;

X0=Points(:);

%only for two cam
X0 = [Camera{1,1};Camera{1,2};Camera{2,1};Camera{2,2};X0];

[X, ~]=HyowonLM(@res_comp,X0,PARAMS);

%options = optimoptions('lsqnonlin','Display','iter','UseParallel',true,...
%                          'Algorithm','levenberg-marquardt',...
%                          'Jacobian', 'on');
                          %'TolX',1e-10,...
                          %'TolFun',1e-10,...
                          %'MaxIter',500);
                      
%X = lsqnonlin(@(X) res_comp(X,PARAMS), X0, [], [], options);

Points = reshape(X(7*PARAMS{2}+1:end),[3 PARAMS{1}]);
Points = Points';
X(1:3*PARAMS{1}) = [];

Camera = zeros(3,4,2);
Camera(1:3,1:3,1) = quat2rotm(X(7:10));
Camera(1:3,1:3,2) = quat2rotm(X(11:14));
Camera(1:3,4,1) = -Camera(1:3,1:3,1)*X(1:3); 
Camera(1:3,4,2) = -Camera(1:3,1:3,2)*X(4:6);



                      
function [F,J] = res_comp(X,PARAMS)
                      
Obs = PARAMS{3};
K   = PARAMS{4};
    
npoints = PARAMS{1}; %number of points
J = zeros(2*PARAMS{2}*npoints,7*PARAMS{2}+3*PARAMS{1});
error = zeros(2, PARAMS{1}*PARAMS{2});

count = 1;
    for i=1:PARAMS{2}
        C = Camera{1,i};
        Q = Camera{2,i};
        R = quat2rotm(Camera{2,i});
        P = K * R * [eye(3) -C];
        
        qw = Q(1);
        qx = Q(2);
        qy = Q(3);
        qz = Q(4);
        
        
        for j=1:size(Obs{i},1)
            
            XIHom = X(7*PARAMS{2} + (1:3) + 3*(Obs{i}(j,3)-1) );
            XHom = [XIHom;1];

            error(count)     = (Obs{i}(j, 1) -  (P(1, :) * XHom) / (P(3, :) * XHom));
            error(count + 1) = (Obs{i}(j, 2) -  (P(2, :) * XHom) / (P(3, :) * XHom));

            
            % Calibration parameters.
            px = K(1, 3);
            py = K(2, 3);
            fx = K(1, 1);
            fy = K(2, 2);
            
            u = [(fx * R(1,1) + px * R(3,1)) (fx * R(1,2) + px * R(3,2)) (fx * R(1,3) + px * R(3,3))] * [XIHom - C];
            v = [(fy * R(2,1) + py * R(3,1)) (fy * R(2,2) + py * R(3,2)) (fy * R(2,3) + py * R(3,3))] * [XIHom - C];
            w = [ R(3,1)                      R(3,2)                      R(3,3)]                     * [XIHom - C];
            
            % partial f / R
            pUpR = [fx * (XIHom - C)'  zeros(1, 3)        px * (XIHom - C)'];
            pVpR = [zeros(1, 3)        fy * (XIHom - C)'  px * (XIHom - C)'];
            pWpR = [zeros(1, 3)        zeros(1, 3)             (XIHom - C)'];
                
            pFpR = [(w * pUpR - u * pWpR) / w^2;
                    (w * pVpR - v * pWpR) / w^2];
                
            % partial q / R (note: quaternion is defined as [qw qx qy qz])
            %{
                   %  w         x        y         z
            pRpQ = [ 0         0        -4 * qy   -4 * qz;
                    -2 * qz    2 * qy    2 * qx   -2 * qw;
                     2 * qy    2 * qz    2 * qw    2 * qx;
                     2 * qz    2 * qy    2 * qx    2 * qw;
                     0        -4 * qx    0        -4 * qz;
                    -2 * qx   -2 * qw    2 * qz    2 * qy;
                    -2 * qy    2 * qz    2 * qw    2 * qx;
                     2 * qx    2 * qw    2 * qz    2 * qy;
                     0        -4 * qx   -4 * qy    0     ];  
            %}
            
            % derivate in order to [qx qy qz qw]
            pRpQ = [ 0        -4 * qy   -4 * qz    0     ;
                     2 * qy    2 * qx   -2 * qw   -2 * qz;
                     2 * qz    2 * qw    2 * qx    2 * qy;
                     2 * qy    2 * qx    2 * qw    2 * qz;
                    -4 * qx    0        -4 * qz    0;
                    -2 * qw    2 * qz    2 * qy   -2 * qx;
                     2 * qz    2 * qw    2 * qx   -2 * qy;
                     2 * qw    2 * qz    2 * qy    2 * qx;
                    -4 * qx   -4 * qy    0         0     ];
            
            % partial f / C
            pUpC = -1 * [(fx * R(1,1) + px * R(3,1)) (fx * R(1,2) + px * R(3,2)) (fx * R(1,3) + px * R(3,3))];
            pVpC = -1 * [(fy * R(2,1) + py * R(3,1)) (fy * R(2,2) + py * R(3,2)) (fy * R(2,3) + py * R(3,3))]; 
            pWpC = -1 * [ R(3,1)                      R(3,2)                      R(3,3)];
            
            pFpC = [(w * pUpC - u * pWpC) / w^2;
                    (w * pVpC - v * pWpC) / w^2];
            
            % partial f / X
            pUpX = [(fx * R(1,1) + px * R(3,1)) (fx * R(1,2) + px * R(3,2)) (fx * R(1,3) + px * R(3,3))];
            pVpX = [(fy * R(2,1) + py * R(3,1)) (fy * R(2,2) + py * R(3,2)) (fy * R(2,3) + py * R(3,3))];
            pWpX = [ R(3,1)                      R(3,2)                      R(3,3)];
            
            pFpX = [(w * pUpX - u * pWpX) / w^2;
                    (w * pVpX - v * pWpX) / w^2];
            
            J(count : (count + 1), :) = [zeros(2, 7 * (i - 1)) pFpR * pRpQ pFpC zeros(2, 7 * (PARAMS{2} - i)) zeros(2, 3 * (j - 1)) pFpX zeros(2, 3*( PARAMS{1} - j ))];

            count = count + 2; 
            
        end
    end
F=error(:);
end

end