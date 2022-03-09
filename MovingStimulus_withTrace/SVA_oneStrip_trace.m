% This code makes cone activation map for the input AOSLO video and shows
% the results on Votonoi diagram (gray scale)
% Asieh Daneshi
% last modification: 05.11.2020
tic
close all
clear
clc
warning('off','all')
%%
% Find_E_January8
FindPRL_November13
clearvars -except alllocs_OneStrip crossFlagOne crossFlagStartOne crossFlagEndOne
% keptFrames=setdiff(1:30,Outliers);   % removing outliers. 30 is the number of frames
keptFrames=find(~isnan(alllocs_OneStrip(:,1)));
PRLX=350;   % PRL location
PRLY=383.5;
clc
%% ========================================================================
fprintf('Please select the cone location file.\n');
[matName,PathName]=uigetfile('BAK*.mat','MultiSelect', 'off');
load([PathName,filesep,matName])
clc
fprintf('Please select the sumnorm image.\n');
[imName,PathName]=uigetfile('BAK*.tiff','MultiSelect', 'off');
Im=imread([PathName,filesep,imName]);
clc
fprintf('Please select the video you aim to analyze. (online stabilized video)\n');
% [videoName,PathName]=uigetfile('*stabilized_del.avi','MultiSelect', 'off');
[videoName,PathName]=uigetfile('*stabilized.avi','MultiSelect', 'off');
myVideo=VideoReader([PathName,filesep,videoName]);
fNom=myVideo.FrameRate*myVideo.Duration;    % Number of frames in the video
ref=im2double(read(myVideo,crossFlagStartOne));
[~,varargout]=corr2d(Im,ref);
[rref,cref]=find(ref);
clc
% Removing conelocs outside boxposition -----------------------------------
conelocs(conelocs(:,1)<boxposition(1),:)=[];
conelocs(conelocs(:,1)>boxposition(1)+boxposition(3),:)=[];
conelocs(conelocs(:,2)<boxposition(2),:)=[];
conelocs(conelocs(:,2)>boxposition(2)+boxposition(4),:)=[];
% Shifting the cone locations to the upper right corner of the image (for
% later use on single frames, not sumnorm) --------------------------------
conelocsN=conelocs;
conelocsN(:,1)=conelocsN(:,1)-boxposition(1)-100;
conelocsN(:,2)=conelocsN(:,2)-boxposition(2)-100;
%% ========================================================================
% =========================================================================
% Initiating variables ----------------------------------------------------
FullField=512;  % In pixel of imaging grid
HalfFF=floor(FullField/2);    % Middle point of FullField
ap_field=FullField;	% The field size, in 'pixels' (typically 512)
% zernike_pupil is the pupil size for the wavefront testing (usually==psf_pupil):
[psf_pupil,zernike_pupil]=deal(7);	% The pupil size for testing
% "deal" simply matches up the input and output lists. For example,
% [a,b,c,...]=deal(5), mean a=5,b=5,c=5
Zoom=10;	% draw PSF much smoother by "zooming in"
% ZF=200;   % additional Zoom Factor (used for cone delivery)
PixPerDegree=600;   % It means each arcmin is 10 pixels
% field_size is the field size, in ARCMIN; (60 arcmin per degree, so 72 for a 1.2deg field):
field_size=(FullField/PixPerDegree*60);    % In 'arcmin'
% field_size=(FullField/PixPerDegree*60)/Zoom;    % In 'arcmin'
lambda=.543;    % The stimulus wavelength;
diff_limited=1; % if diff_limited==1, set all Zernike coefficients to zero; else, input a given .zer from HSWS
defocus=0.03;  % in diopter, 0.01?
% -------------------------------------------------------------------------
% Building PSF
myPSF=GeneratePSF(ap_field,psf_pupil,zernike_pupil,field_size,lambda,diff_limited,defocus);


%% ========================================================================
ZF=2;  % zoom factor
[v,c]=voronoin(conelocs(:,1:2)*ZF);    % N-D Voronoi diagram

A=zeros(length(c),1) ;  % Matrix that will contain area of each voronoi cell

% computing the area of each voronoi cell ---------------------------------
for i=1:length(c)
    v1=v(c{i},1) ;
    v2=v(c{i},2) ;
    A(i)=polyarea(v1,v2);
end
%--------------------------------------------------------------------------
%% ========================================================================
ref=imresize(im2double(read(myVideo,1)),ZF);
sz=size(ref);
S=zeros(sz);

% cropSize is the size of that part of image, around PRL that we are
% focusing at. Please only assign odd numbers to cropSize -----------------
cropSize=250*ZF+1;
conelocIm=zeros(cropSize,cropSize);

croppedFrames=zeros(cropSize,cropSize);



CurrentFrame=imresize(im2double(read(myVideo,1)),ZF);
conelocsNN=[conelocsN(:,1)*ZF+PRLX*ZF-(cropSize-1)/2,conelocsN(:,2)*ZF+PRLY*ZF-(cropSize-1)/2];
% ---------------------------------------------------------------------
croppedFrames=CurrentFrame(PRLX*ZF-(cropSize-1)/2:PRLX*ZF+(cropSize-1)/2,PRLY*ZF-(cropSize-1)/2:PRLY*ZF+(cropSize-1)/2);
clear croppedCL
croppedCL=conelocsNN((conelocsNN(:,1)>=PRLX*ZF-(cropSize-1)/2) & (conelocsNN(:,1)<=PRLX*ZF+(cropSize-1)/2) & (conelocsNN(:,2)>=PRLY*ZF-(cropSize-1)/2) & (conelocsNN(:,2)<=PRLY*ZF+(cropSize-1)/2),:);
croppedConeLocs=croppedCL;
croppedConeLocs(:,1,1)=croppedConeLocs(:,1,1)-PRLX*ZF+(cropSize-1)/2;
croppedConeLocs(:,2,1)=croppedConeLocs(:,2,1)-PRLY*ZF+(cropSize-1)/2;
% croppedConeLocs(:,1,1)=croppedConeLocs(:,1,1)-alllocs_OneStrip(1,1)*ZF;
% croppedConeLocs(:,2,1)=croppedConeLocs(:,2,1)-alllocs_OneStrip(1,2)*ZF;
[croppedVoronoi_v,croppedVoronoi_c]=voronoin(croppedConeLocs(:,1:2,1));    % N-D Voronoi diagram
for n=1:length(croppedVoronoi_c)
    croppedVoronoi_c{n,1}(length(croppedVoronoi_c{n,1})+1)=croppedVoronoi_c{n,1}(1);
end
currentA=zeros(length(croppedVoronoi_c),1);	% Matrix that will contain area of each voronoi cell
for i=1:length(croppedVoronoi_c)
    v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
    v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
    currentA(i)=polyarea(v1,v2);
    if ~isnan(currentA(i))&&(currentA(i)<=250*(ZF)^2)
        currentD(i)=maxDiameter(croppedVoronoi_c{i},croppedVoronoi_v);   % Cone diameter
        SpotSize=floor(2*currentD(i));  % Ask Niklas!
        Aperture=currentD(i)*0.48;   % proportion of inner segement diameter (ISD) that functions a light collection aperture (FWHM) from MacLeod et al, 1992; set range as desired
        cA=Aperture./2.25482;
        HalfSS=floor(SpotSize/2);   % Half of SpotSize
        spot=single(fspecial('gaussian',SpotSize,cA));	%generate Gaussian (SpotSize=ConeSize pixels wide)
        spot=spot./max(spot(:));	%normalize Gaussian
        sz_conelocIm=size(conelocIm);
        % In this loop we go for each cone location in the selected
        % part of the current frame -----------------------------------
        croppedConeLocs1=ceil(croppedConeLocs(i,1,1));
        croppedConeLocs2=ceil(croppedConeLocs(i,2,1));
        if(croppedConeLocs1>0 && croppedConeLocs2>0)    % remove meaningless cone locations
            tempConeIm=zeros(cropSize); % Make a whole black image
            % Whiten the only one pixel in the position of the current cone
            tempConeIm(ceil(croppedConeLocs(i,2,1)),ceil(croppedConeLocs(i,1,1)))=1;
            % Convolve the image of the current cone with a PSF
            filtSpot=conv2(tempConeIm,spot,'same');
            C=zeros(size(filtSpot));
            [rt,ct]=find(filtSpot~=0);
            % Find out which parts of the convolved PSF is inside the current voronoi cell
            inVoronoi=inpolygon(ct,rt,croppedVoronoi_v(croppedVoronoi_c{i},1),croppedVoronoi_v(croppedVoronoi_c{i},2));
            filtSpotVB=zeros(sz_conelocIm(1:2));
            rIn=rt(inVoronoi);
            cIn=ct(inVoronoi);
            for j2=1:length(rIn)
                filtSpotVB(rIn(j2),cIn(j2))=filtSpot(rIn(j2),cIn(j2));
            end
            conelocIm=conelocIm+filtSpotVB;
        end
    end
end

allA(1:length(currentA),1)=currentA;

% preparing an avi file to write videos ===================================
% vWrite=VideoWriter([videoName(1:end-4),'_stim.avi'],'Grayscale AVI');
vWrite=VideoWriter([videoName(1:end-4),'_stim_one.avi']);
vWrite.FrameRate=30;
open(vWrite);
load('Es.mat')      % loading mat file containing all the stimuli
myStim=allE(:,:,str2double(videoName(end-17:end-15)));
% myStim=ones(25);
% myStim(6:10,5:25)=0;
% myStim(16:20,5:25)=0;
myPSF=imresize(myPSF,[cropSize,cropSize]);
a2=0;
% for a1=6
for a1=1:fNom
%     if ismember(a1,keptFrames)
        % for a1=keptFrames(1:10)
        CurrentFrame=imresize(im2double(read(myVideo,a1)),ZF);
        if find(CurrentFrame)
            stimIm=zeros(size(CurrentFrame));
            stimIm(alllocs_OneStrip(a1,2)*ZF-size(myStim,1)/2*ZF:alllocs_OneStrip(a1,2)*ZF+size(myStim,1)/2*ZF,alllocs_OneStrip(a1,1)*ZF-size(myStim,2)/2*ZF:alllocs_OneStrip(a1,1)*ZF+size(myStim,2)/2*ZF)=imresize(myStim,[size(myStim,1)*ZF+1,size(myStim,2)*ZF+1]);  % For square stimulus
            stimIm_translated=stimIm(PRLY*ZF-(cropSize-1)/2:PRLY*ZF+(cropSize-1)/2,PRLX*ZF-(cropSize-1)/2:PRLX*ZF+(cropSize-1)/2);
            if ~isempty(find(stimIm_translated~=0))
                a2=a2+1;
                indStim(a2)=a1;
                % =========================================================================
                % Convolving the stimulus image with the PSF ------------------------------
%                 filtStim=conv2(stimIm_translated,myPSF,'same');
                filtStim=CONVOLVE(stimIm_translated,myPSF);
                % figure;imshow(filtStim,[])
                % Multiplying stimIm and conelocIm(:,:,a2)
                OutputIm=conelocIm.*(filtStim);
%                 OutputImLog=log10(OutputIm);
%                 OutputImLogCropped=zeros(size(OutputImLog));
%                 OutputImLogCropped(400:600,400:600)=OutputImLog(400:600,400:600);
                [rt,ct]=find(filtStim~=0);
                croppedConeLocstemp=croppedConeLocs((croppedConeLocs(:,1)>=min(ct(:))-100) & (croppedConeLocs(:,1)<=max(ct(:))+100) & (croppedConeLocs(:,2)>=min(rt(:))-100) & (croppedConeLocs(:,2)<=max(rt(:))+100),:);
                
                rr=zeros(1,length(croppedConeLocstemp));
                for m=1:length(croppedConeLocstemp)
                    rr(m)=find(croppedConeLocstemp(m,1)==croppedConeLocs(:,1)&croppedConeLocstemp(m,2)==croppedConeLocs(:,2));
                end
                %         croppedConeLocstemp(:,1,1)=croppedConeLocstemp(:,1,1)-min(rt(:));
                %         croppedConeLocstemp(:,2,1)=croppedConeLocstemp(:,2,1)-min(ct(:));
                %         [croppedVoronoi_vtemp,croppedVoronoi_ctemp]=voronoin(croppedConeLocstemp(:,1:2));
                %% ========================================================================
                OutputImF=zeros(cropSize,cropSize);
                %         ii=1;
                for i=rr
                    v1=croppedVoronoi_v(croppedVoronoi_c{i},1);
                    v2=croppedVoronoi_v(croppedVoronoi_c{i},2);
                    currentA(i)=polyarea(v1,v2);
                    if ~isnan(currentA(i))&&(currentA(i)<=250*(ZF)^2)
                        % In this loop we go for each cone location in the selected
                        % part of the current frame -----------------------------------
                        croppedConeLocs1=ceil(croppedConeLocs(i,1));
                        croppedConeLocs2=ceil(croppedConeLocs(i,2));
                        if(croppedConeLocs1>0 && croppedConeLocs2>0)    % remove meaningless cone locations
                            tempConeIm=zeros(cropSize); % Make a whole black image
                            % Whiten the only one pixel in the position of the current cone
                            tempConeIm(ceil(croppedConeLocs(i,1,1)),ceil(croppedConeLocs(i,2,1)))=1;
                            % Convolve the image of the current cone with a PSF
                            %             filtSpot=conv2(tempConeIm,spot,'same');
                            %             C=zeros(size(filtSpot));
                            %                     [rt,ct]=find(OutputIm~=0);
                            
                            
                            % Find out which parts of the convolved PSF is inside the current voronoi cell
                            inVoronoiF=inpolygon(ct,rt,croppedVoronoi_v(croppedVoronoi_c{i},1),croppedVoronoi_v(croppedVoronoi_c{i},2));
                            filtSpotF=zeros(sz_conelocIm(1:2));
                            
                            rIn=rt(inVoronoiF);
                            cIn=ct(inVoronoiF);
                            %                     filtSpotVB(rIn,cIn)=filtSpot(rIn,cIn);
                            for j2=1:length(rIn)
                                filtSpotF(rIn(j2),cIn(j2))=OutputIm(rIn(j2),cIn(j2));
                            end
                            [rz,cz]=find(filtSpotF~=0);
                            %                     if ~isempty(rz)
                            %                         Activity1(a2,ii,1:(max(rz)-min(rz)+1)*(max(cz)-min(cz)+1))=reshape(filtSpotF(min(rz):max(rz),min(cz):max(cz)),[1,(max(rz)-min(rz)+1)*(max(cz)-min(cz)+1)]);
                            %                         ii=ii+1;
                            %                     end
                            mySum=sum(filtSpotF(:));
                            %                     Activity2(a2,i)=mySum;
                            filtSpotVB_sum=zeros(sz_conelocIm(1:2));
                            for j3=1:length(rIn)
                                filtSpotVB_sum(rIn(j3),cIn(j3))=mySum;
                            end
                            OutputImF=OutputImF+filtSpotVB_sum;
                        end
                    end
                end
                %% ========================================================================
                OutImF=OutputImF;
                OutputImF=OutputImF/max(OutImF(:));
                picMat_area(1:cropSize,1:cropSize)=OutputImF;
                picMat_activity(1:cropSize,1:cropSize)=OutputIm;
                save(['picMat_areaOne',num2str(a1),'.mat'],'picMat_area');
                save(['picMat_activityOne',num2str(a1),'.mat'],'picMat_activity');
                %         str1=[num2str(zeros(1,3-numel(num2str(fix(a2))))),num2str(a2)];
                %         str1=str1(find(~isspace(str1)));
                %         str2=[num2str(zeros(1,3-numel(num2str(fix(a2/vWrite.FrameRate))))),num2str(a2/vWrite.FrameRate,'%0.3f')];
                %         str2=str2(find(~isspace(str2)));
                %         OutputImTexted=insertText(OutputImF,[5 2980],...
                %             [str1,':',str2],'FontSize',60,...
                %             'BoxColor','black','BoxOpacity',1,'TextColor','white','AnchorPoint','LeftBottom');
                %         fig1=figure;imshow(OutputImTexted,[])
                fig1=figure;imshow(1-OutputImF,[])
                hold on
                for n=1:length(croppedVoronoi_c)
                    plot(croppedVoronoi_v(croppedVoronoi_c{n},1),croppedVoronoi_v(croppedVoronoi_c{n},2),'-','color','k')
                end
                %             [ra,ca]=find(stimIm_translated~=0);
                %             plot(fix((max(ca(:))-min(ca(:)))/2+min(ca(:))),fix((max(ra(:))-min(ra(:)))/2+min(ra(:))),'r+','markersize',10,'linewidth',2)
                %             plot(conelocsNN(:,1)-PRLX*ZF+1251,conelocsNN(:,2)-PRLY*ZF+1251,'g.')
                plot(nanmedian(alllocs_OneStrip(keptFrames,2))*ZF-PRLY*ZF+1251,nanmedian(alllocs_OneStrip(keptFrames,1))*ZF-PRLX*ZF+1251,'m+','markersize',15,'linewidth',4)
                plot(alllocs_OneStrip(a1,1)*ZF-PRLX*ZF+1251,alllocs_OneStrip(a1,2)*ZF-PRLY*ZF+1251,'r+','markersize',15,'linewidth',3)
                %     text(5,2950,[str1,':',str2],'Color','white','FontSize',10)
                %             str1=[num2str(zeros(1,3-numel(num2str(fix(a1))))),num2str(a1)];
                %             str1=str1(find(~isspace(str1)));
                %             str2=[num2str(zeros(1,4-numel(num2str(fix(a1/vWrite.FrameRate*1000))))),num2str(fix(a1/vWrite.FrameRate*1000))];
                %             str2=str2(find(~isspace(str2)));
                %             annotation('textbox', [.07, .06 0.1 0.02], 'String', [str1,': ',str2],'color','white','fontsize',12,'BackgroundColor','black','FitBoxToText','on')
                myFrame=getframe;
                writeVideo(vWrite,myFrame)
                %             imwrite(myFrame.cdata,['f_',num2str(a1),'.tiff'])
                %             imwrite(1-OutputImF,['frame_',num2str(a1),'.tiff'])
                close(fig1)
            end
%             a2=a2+1;
%         else
%             OutputImF=ones(cropSize,cropSize);
%             fig1=figure;imshow(OutputImF,[])
%             hold on
%             for n=1:length(croppedVoronoi_c)
%                 plot(croppedVoronoi_v(croppedVoronoi_c{n},1),croppedVoronoi_v(croppedVoronoi_c{n},2),'-','color','k')
%             end
%             plot(nanmedian(alllocs_OneStrip(keptFrames,2))*ZF-PRLY*ZF+1251,nanmedian(alllocs_OneStrip(keptFrames,1))*ZF-PRLX*ZF+1251,'m+','markersize',15,'linewidth',4)
%             %     text(5,2950,[str1,':',str2],'Color','white','FontSize',10)
%             str1=[num2str(zeros(1,3-numel(num2str(fix(a1))))),num2str(a1)];
%             str1=str1(find(~isspace(str1)));
%             str2=[num2str(zeros(1,4-numel(num2str(fix(a1/vWrite.FrameRate*1000))))),num2str(fix(a1/vWrite.FrameRate*1000))];
%             str2=str2(find(~isspace(str2)));
%             annotation('textbox', [.07, .06 0.1 0.02], 'String', [str1,': ',str2],'color','white','fontsize',12,'BackgroundColor','black','FitBoxToText','on')
%             myFrame=getframe;
%             writeVideo(vWrite,myFrame)
% %             imwrite(OutputImF,['frame_',num2str(a1),'.tiff'])
%             close(fig1)
%         end
%     else
%         OutputImF=ones(cropSize,cropSize);
%         fig1=figure;imshow(OutputImF,[])
%         hold on
%         for n=1:length(croppedVoronoi_c)
%             plot(croppedVoronoi_v(croppedVoronoi_c{n},1),croppedVoronoi_v(croppedVoronoi_c{n},2),'-','color','k')
%         end
%         plot(nanmedian(alllocs_OneStrip(keptFrames,2))*ZF-PRLY*ZF+1251,nanmedian(alllocs_OneStrip(keptFrames,1))*ZF-PRLX*ZF+1251,'m+','markersize',15,'linewidth',4)
%         %     text(5,2950,[str1,':',str2],'Color','white','FontSize',10)
%         str1=[num2str(zeros(1,3-numel(num2str(fix(a1))))),num2str(a1)];
%         str1=str1(find(~isspace(str1)));
%         str2=[num2str(zeros(1,4-numel(num2str(fix(a1/vWrite.FrameRate*1000))))),num2str(fix(a1/vWrite.FrameRate*1000))];
%         str2=str2(find(~isspace(str2)));
%         annotation('textbox', [.07, .06 0.1 0.02], 'String', [str1,': ',str2],'color','white','fontsize',12,'BackgroundColor','black','FitBoxToText','on')
%         myFrame=getframe;
%         writeVideo(vWrite,myFrame)
% %         imwrite(OutputImF,['frame_',num2str(a1),'.tiff'])
%         close(fig1)
%         %     a2=a2+1;
    end
end
% a2=a2-1;
save('croppedVoronoi_c.mat','croppedVoronoi_c');
save('croppedVoronoi_v.mat','croppedVoronoi_v');
save('keptFramesOne.mat','keptFrames');
save('crossFlagOne.mat','crossFlagOne');
save('alllocs_OneStrip.mat','alllocs_OneStrip');
save('keptFramesOne.mat','keptFrames');
% save('stripShifts1.mat','stripShifts');

% close all
close(vWrite)
%% showing the result with colors on sumNorm image ========================
figure;imshow(zeros(512*ZF,512*ZF),[])
% figure;imshow(imresize(Im,ZF),[])
cmap=flip(parula(100),1);   % Defining a new color map
hold on
for i=1:length(A)
    v1=v(c{i},1);
    v2=v(c{i},2);
    if A(i)<100*(ZF^2)
        patch(v1-boxposition(1,1)*ZF,v2-boxposition(1,2)*ZF,cmap(ceil((floor(A(i))+1)/(ZF^2)),:)); % Fill a voronoi cell with color
    end
end
rectangle('Position',[mean(croppedConeLocs(:,1))+boxposition(1,1)*ZF-cropSize*ZF/2,mean(croppedConeLocs(:,2))+boxposition(1,2)*ZF-cropSize*ZF/2,cropSize*ZF,cropSize*ZF],'EdgeColor','w')
% plot(croppedConeLocs(:,1)+boxposition(1,1)*ZF,croppedConeLocs(:,2)+boxposition(1,2)*ZF,'w.')
% =========================================================================
warning('on','all')
% toc