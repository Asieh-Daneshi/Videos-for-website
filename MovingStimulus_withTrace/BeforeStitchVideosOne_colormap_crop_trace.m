%% This code makes the electrophysiology based videos (videos with traces) 
% from 'picMat_areaOne.mat' files, it must be run after
% 'SVA_oneStrip_trace.m' and before 'StitchVideosOne_colormap_crop_trace.m'
% Asieh Daneshi January. 2021

SVA_oneStrip_trace
clc
ZF=2;
% load('croppedVoronoi_c.mat')
% load('croppedVoronoi_v.mat')
load('crossFlagOne.mat')

xCenter=fix(cellfun(@(index) mean(croppedVoronoi_v(index,1)),croppedVoronoi_c));
yCenter=fix(cellfun(@(index) mean(croppedVoronoi_v(index,2)),croppedVoronoi_c));
centersOld=[xCenter,yCenter];
centers=nan(size(centersOld));
centers((xCenter>0) & (xCenter<=250*ZF+1) & (yCenter>0) & (yCenter<=250*ZF+1),:)=centersOld((xCenter>0) & (xCenter<=250*ZF+1) & (yCenter>0) & (yCenter<=250*ZF+1),:);
fNom=30;
totalActivity=zeros(15,length(centers));
a2=0;

for a1=1:fNom
    if ismember(a1,crossFlagOne)
        a2=a2+1;
        Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
        allPicMat(:,:,a2)=Stimp.picMat_area;
        for a3=1:length(centers)
            if ~isnan(centers(a3,2)) && ~isnan(centers(a3,1))
                picMat_area_centers(1,a3)=Stimp.picMat_area(centers(a3,2),centers(a3,1));
            end
        end
        totalActivity(a1,:)=picMat_area_centers;
    else
        totalActivity(a1,:)=zeros(1,length(centers));
    end
end
totalActivityBounded=zeros(size(totalActivity));
totalActivityBounded(totalActivity>0.02)=totalActivity(totalActivity>0.02);
% figure;imagesc(totalActivity)
% figure;imagesc(totalActivityBounded)

%% ========================================================================
currentActivity=totalActivity;
for a1=3:fNom
    for a2=1:length(totalActivity)
        if totalActivity(a1,a2)==0 && totalActivity(a1-1,a2)==0
            currentActivity(a1,a2)=0;
        else            
            currentActivity(a1,a2)=max(currentActivity(a1,a2),currentActivity(a1-1,a2)/2);
        end
    end
    wholeActivity(a1,:)=currentActivity(a1,:);
end
% Show the results ========================================================
for a=1:30
    % single frames -------------------------------------------------------
    fig1=figure;
%     subplot(1,2,1)
    imshow(ones(250*ZF+1),[])
    hold on
    for i=1:length(centers)
        v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
        v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
        patch(v1,v2,1-totalActivity(a,i),'EdgeColor','none'); % Fill a voronoi cell with color
    end
    myFrame=getframe;
    imwrite(rgb2gray(myFrame.cdata),['SingleFrame_',num2str(a),'.tif']);
    hold off
    close(fig1)
    % Electrophysiology ---------------------------------------------------
% %     figure;
%     subplot(1,2,2)
%     imshow(ones(250*ZF+1),[])
%     hold on
%     for i=1:length(allA)
%         v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
%         v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
%         patch(v1,v2,1-wholeActivity(a,i)); % Fill a voronoi cell with color
%     end
%     hold off
    fig2=figure;
    imshow(ones(250*ZF+1),[])
    hold on
    for i=1:length(centers)
        v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
        v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
        patch(v1,v2,1-wholeActivity(a,i),'EdgeColor','none'); % Fill a voronoi cell with color
    end
    myFrame=getframe;
    imwrite(rgb2gray(myFrame.cdata),['Electrophysiology_',num2str(a),'.tif']);
    asi=im2double(rgb2gray(myFrame.cdata));
    save(['Electrophysiology_',num2str(a),'.mat'],'asi');
    hold off
    close(fig2)
end