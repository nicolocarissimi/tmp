%% Single frame recontrsuction
P = 17;
f = 1;       % frame # you want to reconstruct, try different frames
close all;

load('coods.mat', 'coods'); % file build from from evaluateAP_mAPSin(predidxs)
%coods mpii predictions which have full torso with mapping already done.
% structure index is image id

ind = 21; % selecting image 21
indices = []; % getting people from image 21
for j = 1:length(coods)
    if coods(j).index == ind
        indices = [indices j];
    end
end

%imshow(['C:\Users\pgay\Code\burooj\images\' coods(indices(1)).image])
%ok this image is actually 071324479.jpg
imshow('071324479.jpg')

load skeleton_17Pts
feetInd = [13,17];

figure(3)
hold on;

for k = indices
  xy_MPI = [];  
  xy_MPI(1,:) = coods(k).x;
  xy_MPI(2,:) = coods(k).y;
  xy_MPI = [xy_MPI(1,1:12) 0 xy_MPI(1,13:15) 0;xy_MPI(2,1:12) 0 xy_MPI(2,13:15) 0];
  map = [16 2 3 4 5 6 7 8 1 9 10 11 13 12 14 15 17];% mapping evaluation to 3D paper
  xy_MPI = xy_MPI(:,map);
  xy = xy_MPI;
  xy(:, feetInd) = nan;
  view2D(xy,edges);
end
  



figure(4)  
hold on;
        
for i = indices
    
  xy_MPI = [];  
  xy_MPI(1,:) = coods(i).x;
  xy_MPI(2,:) = coods(i).y;

  xy_MPI = [xy_MPI(1,1:12) 0 xy_MPI(1,13:15) 0;xy_MPI(2,1:12) 0 xy_MPI(2,13:15) 0];


  map = [16 2 3 4 5 6 7 8 1 9 10 11 13 12 14 15 17];
  xy_MPI = xy_MPI(:,map);
  xy = xy_MPI;
  x_torso = xy(1,1);
  y_torso = xy(2,1);
  xy(1,:) = xy(1,:) - xy(1,1);
  xy(2,:) = xy(2,:) - xy(2,1);

%% If feet are missing
  xy(:, feetInd) = nan;
  [cost,~, X2,X2b,R,R1,X3D,r,st] = estimatePoseMissingFeet(xy);   
% figure(1)              
% viewWnS2(xy, X2, X2, edges)

  X3D_R = st*r*R*X3D;
  
  if length(indices) > 1 && i == indices(1)
      ref = X3D_R(2,3);
      r_norm = r; 
  end
  
  if i ~= indices(1)
      ref = ref*r_norm/r;
      shift = ref - X3D_R(2,3);
      X3D_R(:,3) = X3D_R(:,3) + shift;
  end
  X3D_R(1,:) =  X3D_R(1,:) + x_torso;
  X3D_R(2,:) =  X3D_R(2,:) + y_torso;
  r
  view3D(xy,X3D_R,edges)
   %viewWnS2(xy, X3D_R, X3D_R, edges)
end