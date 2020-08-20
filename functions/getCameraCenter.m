function center = getCameraCenter(Pmat)
%get projection matrix and convert it to center of the camera
%P mat 3*4 matrix, [R t]

x =  det([ Pmat(:,2), Pmat(:,3), Pmat(:,4) ]);
y = -det([ Pmat(:,1), Pmat(:,3), Pmat(:,4) ]);
z =  det([ Pmat(:,1), Pmat(:,2), Pmat(:,4) ]);
t = -det([ Pmat(:,1), Pmat(:,2), Pmat(:,3) ]);
center = [ x/t; y/t; z/t ];

end