%% This code stitches each frame of the AOSLO with its corresponding 
% electrophysiology cone map
% Asieh Daneshi January. 2021
warning('off','all')
[videoName,PathName]=uigetfile('*stabilized.avi','MultiSelect', 'off');
myVideo=VideoReader([PathName,filesep,videoName]);
fNom=myVideo.FrameRate*myVideo.Duration;    % Number of frames in the video
%% ========================================================================
% This part of the code crops the input video in a window same as the
% window used in 'SVA_oneStrip' -------------------------------------------
ZF=1;
ZFB=2; % Zoom out to the original size
cropSize=250*ZF+1;
% preparing an avi file to write videos ===================================
vWrite=VideoWriter([videoName(1:end-4),'_Stitch.avi']);
vWrite.FrameRate=30;
open(vWrite);
load('alllocs_OneStrip.mat')
load('croppedVoronoi_c.mat')
load('croppedVoronoi_v.mat')
load('keptFramesOne.mat')
load('crossFlagOne.mat')
Extend=80;
PRLX=345-10;   % PRL location
PRLY=383.5+1;
a2=0;
for a1=1:fNom
    % for a1=6
    if ismember(a1,keptFrames)  % Check if the frame is removed because of its low quality
        CurrentFrame=im2double(read(myVideo,a1));
        CurrentFrameCrop=CurrentFrame(alllocs_OneStrip(a1,2)-(cropSize-1)/2:alllocs_OneStrip(a1,2)+(cropSize-1)/2,alllocs_OneStrip(a1,1)-(cropSize-1)/2:alllocs_OneStrip(a1,1)+(cropSize-1)/2);
        
        %         if exist(['picMat_areaOne',num2str(a1),'.mat'],'file')
        if ismember(a1,crossFlagOne)
            a2=a2+1;
            %                 Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
            Stimp=load(['Electrophysiology_',num2str(a1),'.mat']);
            %                 Stim=imresize((Stimp.picMat_area),1/ZFB);
            Stim=imresize((Stimp.asi),[cropSize,cropSize]);
            [sx,sy]=size(Stim);
            fig1=figure;
            imshow([zeros(cropSize,cropSize),1-Stim(1:sx,1:sy)],[])
%             hold on
%             for n=1:length(croppedVoronoi_c)
%                 plot(croppedVoronoi_v(croppedVoronoi_c{n},1)/ZFB+cropSize+1,croppedVoronoi_v(croppedVoronoi_c{n},2)/ZFB+1,'-','color',[0.5 0.5 0.5],'linewidth',0.1)
%             end
            colormap parula
            myFrame=getframe;
            writeVideo(vWrite,myFrame.cdata)
            close(fig1)
        end
    end
end
close all
close(vWrite)
% fig1=figure;
% imshow(zeros(251),[])
% annotation('textbox', [.4, .7 0.1 0.02], 'String', 'Cone','color','white','fontsize',20,'BackgroundColor','black')
% annotation('textbox', [.35, .55 0.1 0.02], 'String', 'Activation','color','white','fontsize',20,'BackgroundColor','black')
% myFrame=getframe;
% writeVideo(vWrite,myFrame.cdata)
% close(fig1)
%
vWrite=VideoWriter([videoName(1:end-4),'_Stitchz.avi']);
vWrite.FrameRate=10;
open(vWrite);
vRead=VideoReader([videoName(1:end-4),'_Stitch.avi']);

a2=0;
for a1=1:fNom
    % for a1=6
    if ismember(a1,keptFrames)  % Check if the frame is removed because of its low quality
        CurrentFrame=imresize(im2double(read(myVideo,a1)),ZF);
        CurrentFrameCrop=CurrentFrame(PRLY*ZF-(cropSize-Extend-1)/2:PRLY*ZF+(cropSize-Extend-1)/2,PRLX*ZF-(cropSize-Extend-1)/2:PRLX*ZF+(cropSize-Extend-1)/2);
        
        %             if exist(['picMat_areaOne',num2str(a1),'.mat'],'file')
        %                 a2=a2+1;
        %                 Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
        %                 Stim=imresize((Stimp.picMat_area),1/ZFB);
        %                 %             Stim=Stim(1+Extend/2:sx-Extend/2,1+Extend/2:sy-Extend/2);
        
        if ismember(a1,crossFlagOne)
            a2=a2+1;
            %                 Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
            Stimp=load(['Electrophysiology_',num2str(a1),'.mat']);
            %                 Stim=imresize((Stimp.picMat_area),1/ZFB);
            Stim=imresize((Stimp.asi),[cropSize,cropSize]);
            fig1=figure;
            imshow([zeros(251,251),1-Stim],[])
            I=im2double(read(vRead,a2));
            [sx,sy,~]=size(I);
            imshow(I(Extend/2+1:sx-Extend/2,Extend*3/2+1:sy-Extend/2,:));
            hold on
            AOSLO_translated=(CurrentFrameCrop-min(CurrentFrameCrop(:)))/(max(CurrentFrameCrop(:))-min(CurrentFrameCrop(:)));
            imshow(AOSLO_translated,[])
            %             colormap hot
            myFrame=getframe;
            writeVideo(vWrite,myFrame.cdata)
            close(fig1)
        end
    end
end
close all
close(vWrite)







%%
vWrite=VideoWriter([videoName(1:end-4),'_Stitch.avi']);
vWrite.FrameRate=30;
open(vWrite);
load('alllocs_OneStrip.mat')
load('croppedVoronoi_c.mat')
load('croppedVoronoi_v.mat')
load('keptFramesOne.mat')
load('crossFlagOne.mat')
Extend=80;
PRLX=345-10;   % PRL location
PRLY=383.5+1;
a2=0;
for a1=1:fNom
    % for a1=6
    if ismember(a1,keptFrames)  % Check if the frame is removed because of its low quality
        CurrentFrame=im2double(read(myVideo,a1));
        CurrentFrameCrop=CurrentFrame(alllocs_OneStrip(a1,2)-(cropSize-1)/2:alllocs_OneStrip(a1,2)+(cropSize-1)/2,alllocs_OneStrip(a1,1)-(cropSize-1)/2:alllocs_OneStrip(a1,1)+(cropSize-1)/2);
        
        %         if exist(['picMat_areaOne',num2str(a1),'.mat'],'file')
        if ismember(a1,crossFlagOne)
            a2=a2+1;
            %                 Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
            Stimp=load(['Electrophysiology_',num2str(a1),'.mat']);
            %                 Stim=imresize((Stimp.picMat_area),1/ZFB);
            Stim=imresize((Stimp.asi),[cropSize,cropSize]);
            [sx,sy]=size(Stim);
            fig1=figure;
            imshow([zeros(cropSize,cropSize),1-Stim(1:sx,1:sy)],[])
%             hold on
%             for n=1:length(croppedVoronoi_c)
%                 plot(croppedVoronoi_v(croppedVoronoi_c{n},1)/ZFB+cropSize+1,croppedVoronoi_v(croppedVoronoi_c{n},2)/ZFB+1,'-','color',[0.5 0.5 0.5],'linewidth',0.1)
%             end
            colormap gray
            myFrame=getframe;
            writeVideo(vWrite,myFrame.cdata)
            close(fig1)
        end
    end
end
close all
close(vWrite)
% fig1=figure;
% imshow(zeros(251),[])
% annotation('textbox', [.4, .7 0.1 0.02], 'String', 'Cone','color','white','fontsize',20,'BackgroundColor','black')
% annotation('textbox', [.35, .55 0.1 0.02], 'String', 'Activation','color','white','fontsize',20,'BackgroundColor','black')
% myFrame=getframe;
% writeVideo(vWrite,myFrame.cdata)
% close(fig1)
%
vWrite=VideoWriter([videoName(1:end-4),'_Stitchy.avi']);
vWrite.FrameRate=10;
open(vWrite);
vRead=VideoReader([videoName(1:end-4),'_Stitch.avi']);

a2=0;
for a1=1:fNom
    % for a1=6
    if ismember(a1,keptFrames)  % Check if the frame is removed because of its low quality
        CurrentFrame=imresize(im2double(read(myVideo,a1)),ZF);
        CurrentFrameCrop=CurrentFrame(PRLY*ZF-(cropSize-Extend-1)/2:PRLY*ZF+(cropSize-Extend-1)/2,PRLX*ZF-(cropSize-Extend-1)/2:PRLX*ZF+(cropSize-Extend-1)/2);
        
        %             if exist(['picMat_areaOne',num2str(a1),'.mat'],'file')
        %                 a2=a2+1;
        %                 Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
        %                 Stim=imresize((Stimp.picMat_area),1/ZFB);
        %                 %             Stim=Stim(1+Extend/2:sx-Extend/2,1+Extend/2:sy-Extend/2);
        
        if ismember(a1,crossFlagOne)
            a2=a2+1;
            %                 Stimp=load(['picMat_areaOne',num2str(a1),'.mat']);
            Stimp=load(['Electrophysiology_',num2str(a1),'.mat']);
            %                 Stim=imresize((Stimp.picMat_area),1/ZFB);
            Stim=imresize((Stimp.asi),[cropSize,cropSize]);
            fig1=figure;
            imshow([zeros(251,251),1-Stim],[])
            I=im2double(read(vRead,a2));
            [sx,sy,~]=size(I);
            imshow(I(Extend/2+1:sx-Extend/2,Extend*3/2+1:sy-Extend/2,:));
            hold on
            AOSLO_translated=(CurrentFrameCrop-min(CurrentFrameCrop(:)))/(max(CurrentFrameCrop(:))-min(CurrentFrameCrop(:)));
            imshow(AOSLO_translated,[])
            %             colormap hot
            myFrame=getframe;
            writeVideo(vWrite,myFrame.cdata)
            close(fig1)
        end
    end
end
close all
close(vWrite)