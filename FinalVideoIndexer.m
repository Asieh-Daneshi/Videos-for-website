%% This code attaches all the input videos, and makes anootations on them
% Enter frame numbers used in each video (when you run
% 'SVA_oneStrip_trace.m', one of the values it returns is the
% 'crossFlagOne', that is automatically recorded in 'crossFlagOne.mat'. It
% contains the number of frames with cross. Copy all these 'crossFlagOne'
% values for all the target videos in 'frames'
% Also, enter all the PRL positions in 'PRLpositions' ( if you don't aim to
% show PRL, you can comment the corresponding line of the code) and CDC 
% from Jenny's data in CDC
% Asieh Daneshi January. 2021
close all
clear
clc
[videoName,PathName]=uigetfile('*.avi','MultiSelect', 'off');
myVideo=VideoReader([PathName,filesep,videoName]);
fNom=myVideo.FrameRate*myVideo.Duration;

vWrite=VideoWriter([videoName(1:end-4),'_indexed.avi']);
vWrite.FrameRate=30;
open(vWrite);
% frames1=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
frames2=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19];
% PRLpositions=[346,373;309,382;339,374;330,376;343,389.5;343,364;345.5,405.5;330,366.5;316,421.5;346,365;342.5,389;363,423.5;315,349;317,380.5;342,397;340.5,367;344.5,349.5;354,391.5;310.5,370.5;384,369.5];
PRLpositions=[346,373;309,382;339,374;330,376;343,389.5;343,364;345.5,405.5;330,366.5;316,421.5;346,365;342.5,389;363,423.5;315,349;317,380.5;342,397;340.5,367;344.5,349.5;310.5,370.5;384,369.5]; %VLs
% load('BAK1041R_2020_04_20_10_17_53_AOSLO.mat')
% PRLpositions=[378,383.5;329.5,377;350,404;348.5,381;326.5,409.5;338.5,389;356,428.5;362,399;361.5,412;352,396;342.5,381;328,363;374.5,410;310,400.5;335.5,411.5;329.5,387;371.5,387;324,357.5;351,397;360.5,388]; %Es
load('BAK1041R_2020_04_20_09_53_26_AOSLO.mat')
load('BAK1041R1_2020_04_20_09_29_14_AOSLO_788_V006_annotated_JLR_DensityResults_150.mat')
% load('BAK1041R1_2020_04_20_09_29_14_AOSLO_788_V006_annotated_JLR_DensityResults_150.mat')
PRLpositions=[395,383.5];
AllPRLpositions=[395,383.5];
CDC=[3.408732080157712e+02,4.068131188550251e+02];
% for a=1:length(PRLpositions)
%     AllPRLpositions(15*(a-1)+1:15*a,1:2)=repmat(PRLpositions(a,1:2),15,1);
% end
% CDC=nanmean(AllPRLpositions)-PRLpositions(1,:);    
% for b=1:fNom
%     frameNumber(b)=mod(b,15);
%     frameNumber(frameNumber==0)=15;
% end
% frameNumber=frameNumber+4;
for a=1:10
%     st1(a)=fix((a-1)/15)+1
%     shiftX(a)=AllPRLpositions(a,1)-AllPRLpositions(1,1);
%     shiftY(a)=AllPRLpositions(a,2)-AllPRLpositions(1,2);
    I=read(myVideo,a);
%     IText=insertText(I,[50,150],'1 arcmin','FontSize',6,'TextColor','white','AnchorPoint','LeftBottom');
    fig1=figure;
%     imshow(IText,[])
    imshow(I,[])
    hold on
    Xs=370:420;
    Ys=ones(length(Xs),1)*340;
    plot(Xs,Ys,'Color','white','LineWidth',1.5)
    plot(176+(WCentroid20(1,1)-PRLpositions(1,1)),176+(WCentroid20(1,2)-PRLpositions(1,2)),'o','markerfacecolor','w','markeredgecolor','r')
    annotation('textbox', [.515, .20 0.1 0.02], 'String', '5 arcmin','color','white','EdgeColor','none','fontsize',6.5,'fontWeight','bold','FitBoxToText','on')
    str1=[num2str(zeros(1,2-numel(num2str(frames1(fix((a-1)/15)+1))))),num2str(frames1(fix((a-1)/15)+1))];
    str1=str1(find(~isspace(str1)));
    str2=[num2str(zeros(1,2-numel(num2str(frames2(a))))),num2str(frames2(a))];
    str2=str2(find(~isspace(str2)));
    str3=[num2str(zeros(1,3-numel(num2str(fix((frames2(a)-1)/vWrite.FrameRate*1000))))),num2str(fix((frames2(a)-1)/vWrite.FrameRate*1000))];
    str3=str3(find(~isspace(str3)));
    annotation('textbox', [.82 .21 0.1 0.02], 'String', [str1,': ',str2,': ',str3],'color','white','fontsize',6.5,'fontWeight','bold','BackgroundColor','none','EdgeColor','none','FitBoxToText','on')
    str4=[num2str(zeros(1,2-numel(num2str(abs(round(rawData(fix((a-1)/15)+1,2)*6)))))),num2str(abs(round(rawData(fix((a-1)/15)+1,2)*6)))];
    
%     str4=str4(find(~isspace(str4)));
    annotation('textbox', [.80 .89 0.1 0.02], 'String', ['Offset: ',str4,'arcsec'],'color','white','fontsize',6.5,'fontWeight','bold','BackgroundColor','none','EdgeColor','none','FitBoxToText','on')
%     annotation('textbox', [.7 .882 0.1 0.02], 'String', ['Gap: ',str4,' arcsec'],'color','white','fontsize',6.5,'fontWeight','bold','BackgroundColor','none','EdgeColor','none','FitBoxToText','on')
    myFrame=getframe;
    writeVideo(vWrite,myFrame.cdata)
    im=frame2im(myFrame); 
    [imind,cm]=rgb2ind(im,256);
    if a==1 
        imwrite(imind,cm,[videoName(1:end-4),'_indexed.gif'],'gif', 'Loopcount',inf,'DelayTime',1/vWrite.FrameRate); 
    else
        imwrite(imind,cm,[videoName(1:end-4),'_indexed.gif'],'gif','WriteMode','append','DelayTime',1/vWrite.FrameRate); 
    end
    close(fig1)
end



close all
close(vWrite)