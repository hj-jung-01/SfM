function obs = proj(point, Cam, K)


point = K*Cam*[point; 1];
point = point/point(3);

obs = point(1:2);

end

