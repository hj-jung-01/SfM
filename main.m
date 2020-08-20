 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the skeleton code of PA2 in EC5301 Computer Vision.              %
% It will help you to implement the Structure-from-Motion method easily.   %
% Using this skeleton is recommended, but it's not necessary.              %
% You can freely modify it or you can implement your own program.          %
% If you have a question, please send me an email to haegonj@gist.ac.kr    %
%                                                      Prof. Hae-Gon Jeon  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

addpath('functions');
addpath('Jacobian_man_rot');
run(strcat(userpath,'/vlfeat-0.9.21/toolbox/vl_setup'));

%% Define constants and parameters
% Constants ( need to be set )

%RANSAC itr
number_of_iterations_for_5_point    = 1000;

% Thresholds ( need to be set )
threshold_of_distance = 1e-7; 

grow_ransac_itr = 10000;
grow_threshold = 1;

% Matrices
K               = [ 1698.873755 0.000000     971.7497705;
                    0.000000    1698.8796645 647.7488275;
                    0.000000    0.000000     1.000000 ];
                
%% Feature extraction and matching
% Load images and extract features and find correspondences. Fill num_Feature, 
% Feature, Descriptor, num_Match and Match hints : use vl_sift to extract features 
% and get the descriptors. use vl_ubcmatch to find corresponding matches between 
% two feature sets.


                
imgPath = 'data_30';
imgType = '/*.jpg';

images_dir  = dir([imgPath imgType]);
    
N = length(images_dir);
cel = cell(N,1);

    for idx = 1:N
        cel(idx,1) =  cellstr(images_dir(idx).name);
        
    end
    
    cel =natsortfiles(cel);
   
    for idx = 1:N
        images_dir(idx).name = char(cel(idx));
    end
    
        images{N,1} = [];

    for idx = 1:N
        images{idx} = imread(fullfile(imgPath,images_dir(idx).name));
    end


D{N,1} = []; % Descriptor Cell 
F{N,1} = []; % Feature Cell
MF{N,N} = []; % Best Match Cell
M = zeros(N,N);
Scores{N,N} = [];

peak_thres = 3;
edge_thres = 10;
foctave = 0;

for i = 1 : N
    im = single(rgb2gray(images{i}));
    %[f,d] = vl_sift(im,'PeakThresh',peak_thres,'EdgeThresh',edge_thres,'FirstOctave',foctave);
    [f,d] = vl_sift(im);
    F{i}=f;        
    D{i}=d;
end

for i = 1 : N
    for j = 1 : N
        if i~= j
        [MF{i,j}, Scores{i,j}] = vl_ubcmatch(D{i}, D{j});
        end
     end
end

    for i = 1 : N
        for j = 1 : N
            M(i,j) =length(MF{i,j}) ;
            if i==j
                M(i,j) = 0;
            end
         end
    end
    
    %Save some features
    %%save('Features','F','D','M','MF','images')


                    
%%

sumc = sum(M');
[argvalue, argmax] = max(sumc(:));
[a , b] =  max(M(argmax,:));
image_sort{1} = [argmax b];

old = image_sort{1};

for image_step=3:length(images)

feats = zeros(1,length(old));
maxfeats = 0;
for j = 1 : length(M)
    if ~ismember(j,old)
        for i = 1: length(old)
            feats(i) = M(j,old(i));
        end
        featsS = sum(feats);
        if maxfeats<featsS
            maxIm = j;
            [arg,ind]= max(feats);
            pairIm = old(ind);
        end
        maxfeats = max(maxfeats,featsS);
    end
end

old = [old maxIm];
image_sort{image_step-1} = [maxIm,pairIm];

end

%% Initialization step
% Estimate E using 8,7-point algorithm or calibrated 5-point algorithm and RANSAC


    
pair = image_sort{1};

ref_feat = MF{pair(1),pair(2)}; % All features of pair 1 and pair 2

% x and y co-ord of all the normalized features in image 1 and imag2 
f_norm1 = inv(K)*[F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
f_norm2 = inv(K)*[F{pair(2)}(1:2,:); ones(1, size(F{pair(2)}(1:2,:),2))];


% x and y cord of matching features  
Q1_feat = f_norm1(:,ref_feat(1,:));
Q2_feat = f_norm2(:,ref_feat(2,:));


Q_feat = [Q1_feat;Q2_feat];
ninliers = 0;     % Number of inliers

% Starting code for RANSAC and Calibrated five point
    for i = 1: number_of_iterations_for_5_point

        % Sampling to get 5 points repeatedly
        ind = randsample(size(Q_feat,2), 5);

        Q1 = Q_feat(1:3,ind);
        Q2 = Q_feat(4:6,ind);

        Evec = calibrated_fivepoint(Q1,Q2);

        nsol = size(Evec,2);
        Emat = permute(reshape(Evec, 3, 3, nsol), [2,1,3]); % Reshaping to get each solution
        E = mat2cell(Emat, 3, 3, ones(1, nsol));

        x1 = Q1_feat;    % Extract x1 and x2 from x
        x2 = Q2_feat;

        nE = length(E);   % Number of solutions to test
        
        for k = 1:nE
          Ex1 = E{k}*x1;
          Etx2 = E{k}'*x2;     

          x2tEx1 = sum(x2.*Ex1);

          d =  x2tEx1.^2 ./ ...
              (Ex1(1,:).^2 + Ex1(2,:).^2 + Etx2(1,:).^2 + Etx2(2,:).^2);

          inliers = find(abs(d) < threshold_of_distance);     % Checking all inliers

          if length(inliers) > ninliers   % Storing the best solution till now
            ninliers = length(inliers);
            bestE = E{k};
            bestInliers = inliers;
          end
        end
    end
    
% Perform Singular value decomposition
[u,s,v]=svd(bestE);
diag_110 = [1 0 0; 0 1 0; 0 0 0];
newE = u*diag_110*v';
[u,s,v] = svd(newE);
w = [0 -1 0; 
     1  0 0;
     0  0 1];
 
% Find Four Possible Solutions  
R1 = u*w*v';
R2 = u*w'*v';
t1 = u(:,3);%./max(abs(u(:,3)));%./max(abs(u(:,3)));
t2 = -u(:,3);%./max(abs(u(:,3)));%./max(abs(u(:,3)));

P1 =   [R1 t1];
P2 =   [R1 t2];

P3 =   [R2 t1];
P4 =   [R2 t2];

P = {P1,P2,P3,P4};

Q1 = Q_feat(1:3,:);
Q2 = Q_feat(4:6,:);
CP1 = eye(3,4);
CP2 = P1;

Xtemp=[];

% Perform the traingulation step to show that not all results are correct
for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X41 =  Xtemp;
X41 = X41(1:3,:);


CP2 = P2;
for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X42 =  Xtemp;
X42 = X42(1:3,:);


CP2 = P3;
for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
end

Xtemp = Xtemp./Xtemp(4,:);
X43 =  Xtemp;
X43 = X43(1:3,:);

CP2 = P4;

for i=1:size(Q_feat,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
    
end

Xtemp = Xtemp./Xtemp(4,:);
X44 =  Xtemp;
X44 = X44(1:3,:);

X_exist = {X41,X42,X43,X44};

P ={P1,P2,P3,P4};
CP1 = [eye(3,3),zeros(3,1)];

Depth = [];

% Now Checking which P is correct by finding camera depth for both cameras

for j = 1:4
X = X_exist{1,j}; 

CP2 = P{1,j};

x =  det([ CP1(:,2), CP1(:,3), CP1(:,4) ]);
y = -det([ CP1(:,1), CP1(:,3), CP1(:,4) ]);
z =  det([ CP1(:,1), CP1(:,2), CP1(:,4) ]);
t = -det([ CP1(:,1), CP1(:,2), CP1(:,3) ]);
c1 = [ x/t; y/t; z/t ]; % finding camera center

CP2 = P{1,j};

x =  det([ CP2(:,2), CP2(:,3), CP2(:,4) ]);
y = -det([ CP2(:,1), CP2(:,3), CP2(:,4) ]);
z =  det([ CP2(:,1), CP2(:,2), CP2(:,4) ]);
t = -det([ CP2(:,1), CP2(:,2), CP2(:,3) ]);
C2{1,j} = [ x/t; y/t; z/t ];  % finding camera center

% Now checking Which cameras have positive depth 
rot2 = CP2(:,1:3); %taking only the rotation matrix
rot1 = CP1(:,1:3);
c2 = C2{1,j};

    for i=1:size(Q_feat,2)
   % Camera 1
    wc1 = rot1(3,:) * (X(1:3,i) - c1(1:3,:));
    depthc1 = (sign(det(rot1)) * wc1) / 1 * norm(rot1(3,:));
    DC1(i) = depthc1(1,1);
    
   %Camera 2
    wc2 = rot2(3,:) * (X(1:3,i) - c2(1:3,:));
    depthc2 = (sign(det(rot2)) * wc2) / 1 * norm(rot2(3,:)); 
    DC2(i) = depthc2(1,1);
    
    end
    D1 = sum(sign(DC1));
    D2 = sum(sign(DC2));
    
    Depth = [Depth D1+D2];


end
[arg,val] = max(Depth);

X = cell2mat(X_exist(val));
X = [X;ones(1,length(X))];
Cmat = cell2mat(P(val));

Pcell{1,length(MF)} = [];
Pcell{1,pair(1)} = [eye(3,3) zeros(3,1)]; 
Pcell{1,pair(2)} = Cmat;

X = X(1:3,:);


% Extracting feature one more time for color extraction
f_unnorm1 = [F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
Q_unnorm = round(f_unnorm1(:,ref_feat(1,:)));
colornew = zeros(3,length(Q_unnorm)); 

for i = 1 : length(Q_unnorm)
    color = images{pair(1)};%(Q_unnorm(2,:),Q_unnorm(1,:),:);
    colornew(:,i) = [color(Q_unnorm(2,i),Q_unnorm(1,i),1), color(Q_unnorm(2,i),Q_unnorm(1,i),2) ,color(Q_unnorm(2,i),Q_unnorm(1,i),3) ];
end


X_exist = X(:,find(X(3,:)>=0));
Q_exist = Q_feat(:,find(X(3,:)>=0));
ref_exist = ref_feat(:,find(X(3,:)>=0));
color_exist = colornew(:,find(X(3,:)>=0));


x2dX3d{1,length(MF)} = [];
x2dX3d{1,pair(1)} =[Q_exist;X_exist;ref_exist(1,:);color_exist]; 
x2dX3d{1,pair(2)} =[Q_exist;X_exist;ref_exist(2,:);color_exist];

X_final = [X_exist; color_exist/255];
SavePLY('2_views.ply', X_final);

%%
%growing step

for growing_step=2:length(image_sort)
    
    pair = image_sort{growing_step};
    pair = [pair(2) pair(1)]; %old new
    
    ref_feat = MF{pair(1),pair(2)};
    
    
    f_norm1 = inv(K)*[F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
    f_norm2 = inv(K)*[F{pair(2)}(1:2,:); ones(1, size(F{pair(2)}(1:2,:),2))];
    
    
    
    old_feat = x2dX3d{1,pair(1)}(10,:);
    [tf1, idx1] = ismembertol(old_feat,ref_feat(1,:));
    ind2 = find(tf1(1,:));
    ref_feat_match = old_feat(1,ind2);
        
    xd2_check = f_norm1(:,ref_feat_match(1,:)); 
    ref_feat_match1 = ref_feat(2,idx1(idx1>0));
    xd2 = f_norm2(:,ref_feat_match1(1,:));
    Xd3 = x2dX3d{1,pair(1)}(7:9,ind2);
    
    Xd3 = [Xd3;ones(1,size(Xd3,2))];
    if(size(Xd3,2)<3)
        break;
    end
    xd2Xd3 = [xd2;Xd3]';
    
    n_grow_inliers = 0;
    
    for itr=1:grow_ransac_itr
        
        
        three_point_ind = randsample(size(xd2Xd3,1), 3);
        
        data_p3p = xd2Xd3(three_point_ind,:);
        
        
        RT=[];
        RT=PerspectiveThreePoint(data_p3p(:,1:6));
        
        nRT=size(RT,1)/4;
        

        for j=1:nRT
            Pcheck = RT(4*(j-1)+1:4*(j-1)+3,:);
            
            KP1X = K*Pcheck*Xd3;
            KP1X = KP1X./KP1X(3,:);
                
            KP2X = K*(Pcell{1,pair(1)})*Xd3;
            KP2X = KP2X./KP2X(3,:);
            
            temp1 = (K*xd2 - KP1X).^2;
            temp2 = (K*xd2_check - KP2X).^2;
            
            temp1 = temp1(1,:)+temp1(2,:);
            temp2 = temp2(1,:)+temp2(2,:);

            dsq = temp1 +temp2;
            d_grow = sqrt(dsq);
            
            grow_inliers = find(abs(d_grow) < grow_threshold);
            
            if length(grow_inliers) > n_grow_inliers
                n_grow_inliers = length(grow_inliers);
                bestP = Pcheck;
                best_grow_inliers = grow_inliers;
            end
        end
                 
    end
    
    Pcell{1,pair(2)} = bestP;
    
    Q1_feat = f_norm1(:,ref_feat(1,:));
    Q2_feat = f_norm2(:,ref_feat(2,:));

    
    Xtemp=[];
    Q1 = [];
    Q2 = [];

    CP1 = Pcell{1,pair(1)};
    CP2 = Pcell{1,pair(2)};
    Q1 = Q1_feat;
    Q2 = Q2_feat;
    
    for i=1:size(Q1,2)
        A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
        [U,S,V] = svd(A);
        Xtemp(:,i) = V(:,4);
    end

    Xtemp = Xtemp./Xtemp(4,:);
    X_grow = Xtemp;
    X_grow = X_grow(1:3,:);
    
    f_unnorm1 = [];
    Q_unnorm = [];
    colornew = [];
    
    f_unnorm1 = [F{pair(1)}(1:2,:); ones(1, size(F{pair(1)}(1:2,:),2))];
    Q_unnorm = round(f_unnorm1(:,ref_feat(1,:)));
    colornew = zeros(3,length(Q_unnorm)); 

    for i = 1 : length(Q_unnorm)
        color = images{pair(1)};%(Q_unnorm(2,:),Q_unnorm(1,:),:);
        colornew(:,i) = [color(Q_unnorm(2,i),Q_unnorm(1,i),1), color(Q_unnorm(2,i),Q_unnorm(1,i),2) ,color(Q_unnorm(2,i),Q_unnorm(1,i),3) ];
    end
    
    X_exist = [];
    Q_exist = [];
    ref_exist = [];
    color_exist = [];
    
    X_exist = X_grow(:,find(X_grow(3,:)>=0));
    Q_exist = [Q1_feat(:,find(X_grow(3,:)>=0)); Q2_feat(:,find(X_grow(3,:)>=0))];
    ref_exist = ref_feat(:,find(X_grow(3,:)>=0));
    color_exist = colornew(:,find(X_grow(3,:)>=0));
    
    x2dX3d{1,pair(2)} =[Q_exist;X_exist;ref_exist(2,:);color_exist];
    
    test1 = round(X_exist(1:3,:)*5);
    [a,indx,c]= unique(test1','rows');
    X_exist = X_exist(1:3,indx);
    color_exist = color_exist(:,indx);
    

    X_add = [];
    X_add = [X_exist; color_exist/255];
    X_final = [X_final X_add];
    
end
Xp1 = X_final(1:3,:);
rot = [-1  0  0;  0 -1  0; 0  0  1];
Xp1 = rot * Xp1;
Xp1 = Xp1(:,find(Xp1(3,:)>=0));
colorp1 = X_final(4:6,:);
colorp1 = colorp1(:,find(Xp1(3,:)>=0));

% Save 3D points to PLY
SavePLY('32_views.ply', X_final);

%%
%bundle adjustment

X_final = [X_exist; color_exist/255];
SavePLY('before_BA.ply', X_final);
%{

%pair = [7,6];
%bundle adjustment for 2 views
%for i=1:2
%    Camera(:,:,i) = [Pcell{pair(i)}; 0 0 0 1];
    %Obs{i} = [x2dX3d{1,pair(2)}(1+3*(i-1):2+3*(i-1),:); (1:length(x2dX3d{1,pair(2)}))]';
%    Q_temp = K*Q_exist(3*i-2:3*i,:);
%    Obs{i} = [Q_temp(1:2,:); (1:length(x2dX3d{1,pair(2)}))];
%end
%Points = x2dX3d{pair(2)}(7:9,:)';
%color_exist = x2dX3d{pair(2)}(11:13,:);

%beforeBA = [Points'; color_exist/255];

%%%only for opt=3, quat
pair = [7,6];
for i=1:2
    Camera{1,i} = getCameraCenter(Pcell{pair(i)});
    Camera{2,i} = rotm2quat(Pcell{pair(i)}(:,1:3))';
end
Points = X_exist;
visible = (1:length(Q_exist));
Obs{1} = vertcat(Q_exist(1:2,:), visible)';
Obs{2} = vertcat(Q_exist(4:5,:), visible)';

opt = 3;
[Points, Camera] = bundleAdj(Points, Camera, Obs, K, opt);

%}


Q_feat1 = K*Q_exist(1:3,:);
Q_feat2 = K*Q_exist(4:6,:);
feats = [Q_feat1(1:2,:); Q_feat2(1:2,:)];
dist = [0 0 0 0 0];

paramters.h = F_vdist(feats, K, dist, [1296,1944,3]);
paramters.K = K; 
paramters.feats = feats;
paramters.rol = 0;
paramters.max_iter = 100;

cam = [rotm2eul(Pcell{pair(2)}(:,1:3))'; Pcell{pair(2)}(:,4);];
x0 = [cam; X_exist(:)];
[x, Residual, J_LM]=HyowonLM(@Bundlecost,x0,paramters);

Points = reshape(x(7:end), 3, [])';
Camera = zeros(3,4,2);
Camera(1:3,:,1) = Pcell{pair(1)};
Camera(1:3,1:3,2) = eul2rotm(x(1:3));
Camera(1:3,4,2) = x(4:6);

bundle_X_points_only = [Points'; color_exist/255];
            
%camera
Xtemp=[];
CP1 = Camera(1:3,:,1);
CP2 = Camera(1:3,:,2);
Q1 = Q_exist(1:3,:);
Q2 = Q_exist(4:6,:);
for i=1:size(Q1,2)
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
    
end

Xtemp = Xtemp./Xtemp(4,:);
bundle_X_BC =  Xtemp;
bundle_X_BC = bundle_X_BC(1:3,:);

bundle_X_camera_only = [bundle_X_BC; color_exist/255];

%camera and point
Xtemp=[];
for i=1:size(Q1,2)
    Q1(:,i)=[proj(Points(i,1:3)', CP1, K);1];
    Q2(:,i)=[proj(Points(i,1:3)', CP2, K);1];
    A = [ Q1(1,i)*CP1(3,:) - CP1(1,:); Q1(2,i)*CP1(3,:) - CP1(2,:); Q2(1,i)*CP2(3,:) - CP2(1,:);  Q2(2,i)*CP2(3,:) - CP2(2,:) ];
    [U,S,V] = svd(A);
    Xtemp(:,i) = V(:,4);
    
end

Xtemp = Xtemp./Xtemp(4,:);
bundle_X_BCP =  Xtemp;
bundle_X_BCP = bundle_X_BCP(1:3,:);

bundle_X_camera_and_pt = [bundle_X_BCP; color_exist/255];

SavePLY('2_views_BP.ply', bundle_X_points_only);
SavePLY('2_views_BC.ply', bundle_X_camera_only);
SavePLY('2_views_BCP.ply', bundle_X_camera_and_pt);




addpath('harris');
imgPath = 'data_30';
imgType = '/*.jpg';

images_dir  = dir([imgPath imgType]);
    
N = length(images_dir);
cel = cell(N,1);

    for idx = 1:N
        cel(idx,1) =  cellstr(images_dir(idx).name);
        
    end
    
    cel =natsortfiles(cel);
   
    for idx = 1:N
        images_dir(idx).name = char(cel(idx));
    end
    
        images{N,1} = [];

    for idx = 1:N
        images{idx} = imread(fullfile(imgPath,images_dir(idx).name));
    end
 
    pair = [11,12]
    
            [~,v,u] = harris(rgb2gray(images{pair(1)}),2,2,3,0);
            feats = [u,v]';
            Nf = size(feats, 2);
            Nf_prev = Nf;

            tracker = vision.PointTracker('MaxBidirectionalError', 0.1);
            initialize(tracker, feats', im2double(images{pair(1)}));
            valid = ones(Nf, 1);
            Ni=length(pair);
            for ni = 2:Ni
               ImgT = im2double(images{pair(ni)});
               [feats_i, valid_i] = step(tracker, ImgT);
               feats(1+(ni-1)*2 : 2+(ni-1)*2, :) = feats_i';
               valid = valid + valid_i;
            end

            falsevalid = valid < Ni;
            feats(:, falsevalid) = [];
            Nf = size(feats,2);
            
            dist = [0 0 0 0 0];
paramters.h = F_vdist(feats, K, dist, [1296,1944,3]);
paramters.K = K; 
paramters.feats = feats;
paramters.rol = 0;
paramters.max_iter = 100;

            X = K\[feats(1,:); feats(2,:); ones(1,Nf)]*100;    
            x0 = [zeros(6*(Ni-1),1); X(:)];

[x, Residual, J_LM]=HyowonLM(@Bundlecost,x0,paramters);
%OPTIONS = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',  'Display', 'iter','Jacobian', 'on');
%[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@(x)Bundlecost(x,paramters),x0,[],[],OPTIONS); % Invoke optimizer

npoints = length(X);
X_ba = reshape(x(7:6+3*npoints),[3 npoints]);
color_exist = zeros(3,npoints).*255;

X_final = [X_ba; color_exist/255];
SavePLY('afterBA.ply', X_final);


function h = F_vdist(feats, K, distort, Imgsize)

    Nf=size(feats,2);
    Ni=size(feats,1)/2;
    u_feat=feats(1:2:end,:);
    v_feat=feats(2:2:end,:);

    k1=distort(1); k2=distort(2); k3=distort(3);
    p1=distort(4); p2=distort(5);

    uv1_undist=[u_feat(:)';v_feat(:)';ones(1,Ni*Nf)];
    nxny1=K\uv1_undist;
    nx_undist=nxny1(1,:);
    ny_undist=nxny1(2,:);
    r_undist=sqrt(nx_undist.^2+ny_undist.^2);
    nx_dist=nx_undist.*(1 + k1*r_undist.^2 + k2*r_undist.^4 + k3*r_undist.^6) + (2*p1*nx_undist.*ny_undist + p2*(r_undist.^2 + 2*nx_undist.^2));
    ny_dist=ny_undist.*(1 + k1*r_undist.^2 + k2*r_undist.^4 + k3*r_undist.^6) + (p1*(r_undist.^2 + 2*ny_undist.^2) + 2*p2*nx_undist.*ny_undist);

    nxny1_dist=[nx_dist;ny_dist;ones(1,Ni*Nf)];
    uv1_dist=K*nxny1_dist;
    u_dist=reshape(uv1_dist(1,:),size(u_feat));
    v_dist=reshape(uv1_dist(2,:),size(v_feat));
    h= ((v_dist-1).*Imgsize(2)+u_dist-1)/(Imgsize(1)*Imgsize(2));
    
end