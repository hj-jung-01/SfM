function [Points, Camera] = bundleAdj(Points, Camera, Obs, K, opt)

if opt == 1
    [Points, Camera] = bundleAdj_lsq(Points, Camera, Obs, K);
elseif opt ==2
    [Points, Camera] = bundleAdj_rot(Points, Camera, Obs, K);
elseif opt ==3
    [Points, Camera] = bundleAdj_quat(Points, Camera, Obs, K);
end

end