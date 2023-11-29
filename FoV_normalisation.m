function IM_mat = FoV_normalisation(IM_mat)

% This function nomorlizes FoV of images with different camera focal plane
% It should be invoked as:
%
%   IM_mat = FoV_normalisation(IM_mat)
%
% where 
%   IM_mat,  is image group from cameras with different focal plane
%           
%   Method, is using bolb detection to estimate scale and descriptor to
%           estiamte corresponding points and directions
%
%  Zhenzhe Han
%  12/2022

len = size(IM_mat,3);
for i = 1:len-1


    Fov_1 = IM_mat(:,:,i);
    Fov_2 = IM_mat(:,:,i+1);

%     Fov_1_e = edge(imgaussfilt(Fov_1,5),'canny');
%     Fov_2_e = edge(imgaussfilt(Fov_2,5),'canny');
    
%     ttl = sprintf( ['\\fontsize{10}',...
%      '{\\color{green}%s} original         ',...    
%      '{\\color{magenta}%s} distorted'],... 
%      char(9632), char(9632));
% 
%  
%     figure('WindowStyle','docked'), imshowpair(Fov_1_e, Fov_2_e)
%     
%     title(ttl);

    ptsFov_1 = detectSURFFeatures(Fov_1,"MetricThreshold",1e-5);
    ptsFov_2 = detectSURFFeatures(Fov_2,"MetricThreshold",1e-5);

    [featuresFov_1,validPtsFov_1] = extractFeatures(Fov_1,ptsFov_1);
    [featuresFov_2,validPtsFov_2] = extractFeatures(Fov_2,ptsFov_2);

    indexPairs = matchFeatures(featuresFov_1,featuresFov_2);

    matchedFov_1 = validPtsFov_1(indexPairs(:,1));
    matchedFov_2 = validPtsFov_2(indexPairs(:,2));
% 
%     figure('WindowStyle','docked')
%     showMatchedFeatures(Fov_1 ,Fov_2 ,matchedFov_1 ,matchedFov_2);
%     title('Putatively matched points (including outliers)');

    [tform, inlierIdx] = estgeotform2d(matchedFov_2,matchedFov_1 ,'similarity');
    inlierFov_2 = matchedFov_2(inlierIdx,:);
    inlierFov_1 = matchedFov_1(inlierIdx,:);
    
%     figure('WindowStyle','docked')
%     showMatchedFeatures(Fov_1,Fov_2,inlierFov_1,inlierFov_2);
%     title('Matching points (inliers only)');
%     legend('ptsOriginal','ptsDistorted');

    invTform = invert(tform);
    Ainv = invTform.A;
    
    ss = Ainv(1,2);
    sc = Ainv(1,1);
    scaleRecovered = hypot(ss,sc);
    disp(['Recovered scale: ', num2str(scaleRecovered)])
    
    % Recover the rotation in which a positive value represents a rotation in
    % the clockwise direction.
    thetaRecovered = atan2d(-ss,sc);
    disp(['Recovered theta: ', num2str(thetaRecovered)])
    
%     outputView = imref2d(size(Fov_1_e));
%     Fov_2_e = imwarp(Fov_2_e,tform,OutputView=outputView);

%     figure('WindowStyle','docked'), imshowpair(Fov_1_e,Fov_2_e)
% 
%     ttl = sprintf( ['\\fontsize{10}',...
%          '{\\color{green}%s} original         ',...    
%          '{\\color{magenta}%s} recovered'],... 
%          char(9632), char(9632));
%      title(ttl);



    outputView = imref2d(size(Fov_1));
    recovered = imwarp(Fov_2,tform,OutputView=outputView);
%     lap=[-1 -1 -1; -1 8 -1; -1 -1 -1];
%     recovered = recovered+1*imfilter(recovered,lap,'replicate');
    IM_mat(:,:,i+1) = recovered;
    
end

end