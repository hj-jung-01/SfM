% Author: Riccardo Giubilato, Padova (Italy) 2016
% mail:   riccardo.giubilato@gmail.com
% https://www.researchgate.net/profile/Riccardo_Giubilato
% -------------------------------------------------------
% Written and tested on Matlab 2016a

function [Points, Camera] = bundleAdj(Points, Camera, Obs, K)
% Using lsqnonlin Matlab function.
% Data structure:
% X=| x y z     x y z     x y z   ..   eul1 eul2 eul3 tx ty tz ... |
%    <point1>  <point2>  <point3> .. < camera ------------ > 
PARAMS{1}=size(Points,1);
PARAMS{2}=size(Camera,3);
PARAMS{3}=Obs;
PARAMS{4}=K;


X0=reshape(Points',1,3*PARAMS{1});

Xcam = []
for i=1:size(Camera,3)
   % Fill with i-th camera parameters
   Xcam = [Xcam rotm2eul(Camera(1:3,1:3,i)) Camera(1:3,4,i)'];

end

X0 = [Xcam X0];
    

options = optimoptions('lsqnonlin','Display','iter','UseParallel',true,...
                          'Algorithm','levenberg-marquardt',...
                          'Jacobian', 'on',...
                          'MaxIter',10);
                      %,...
                      %    'ScaleProblem', 'jacobian',...
                      %    'Jacobian', 'on'
                          %'TolX',1e-10,...
                          %'TolFun',1e-10,...
                          %'MaxIter',500);
                       

X = lsqnonlin(@(X) res_comp(X,PARAMS), X0, [], [], options);

% Ricopio nuove soluzioni
Points = reshape(X(1:3*PARAMS{1}),[3 PARAMS{1}]);
Points = Points';
X(1:3*PARAMS{1}) = [];
for i=1:size(Camera,3)
   Camera(1:3,1:3,i) = eul2rotm(X(1:3)); 
   Camera(1:3,4,i)   = X(4:6)';
   X(1:6)            = [];
end

end

% Calcolo residui
function [F,J] = res_comp(X,PARAMS)

F = [];
Obs = PARAMS{3};
K   = PARAMS{4};

npoints = length(Obs{1}); %number of points
J   = zeros(2*PARAMS{2}*npoints, 6*PARAMS{2}+3*npoints);

for i=1:PARAMS{2}
    for j=1:size(Obs{i},1)
                            % (u,v) of j-th observation in frame i-th
        point = X((1:3) + 3*(Obs{i}(j,3)-1) + 6*PARAMS{2})';
        rotat = eul2rotm(X((1:3)+6*(i-1)));
        trans = X((4:6)+6*(i-1))';
        
        %camparam = [K(1,1), K(2,2), K(1,3), K(2,3), rotm2eul(rotat), trans'];
        camparam = [rotm2eul(rotat), trans'];
        
        F = [F reproj(Obs{i}(j,1:2), point, [rotat trans], K)];
        
        %J_A = computeJacobianCameraParam(point',camparam);
        %J_B = computeJacobianPoint(point',camparam);
        J_A = JacobianCam(point',camparam, K);
        J_B = JacobianPoint(point',camparam,K);
            
        J(2*(j-1)+1+2*npoints*(i-1):2*(j-1)+2+2*npoints*(i-1), 1+6*(i-1):6+6*(i-1)) = J_A;
        J(2*(j-1)+1+2*npoints*(i-1):2*(j-1)+2+2*npoints*(i-1), 12+3*(j-1)+1:12+3*(j-1)+3) = J_B;        
        
    end
end

end
