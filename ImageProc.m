function [m2,intensity ] = ImageProc( filename, thres,sample, frame, fishnum)
%IMAGEPROC default program to process fluorescence imagess
%   sample is index for samples
%   frame is the index for image of each sample
%   fishnum is the number of the image
%   thres is the threshold model to get cell masks
%   2017-3-15
%   the next question is how to subtract background without affect the spot
%   intensity.

%   image processing is divided into two steps, in the first step,
%   background fluorescence is subtracted from each cell. In the second
%   step, the spot instensity is collected

%   para
% mask from fluorescence image
% spot area threshold
% spot intensity threshold
% parameter (0.006, 10, 50, 400) for Organoid2010_c0.tif (turbo 0IPTG),
% thres = 0.015
% parameter (0.006, 10, 1000, 200) for turbo 935
% 2017-10-29
% parameter (0.006, 10,200,200) for dh10B 935
% parameter (0.004, 10, 50, 400) for Organoid1984_c0.tif (DH10B 0IPTG),
% thres = 0.008

% 2017-11-05 parameter (0.07,10,50,1000) for 1105data
% 2017-11-06 parameter (0.019, 10, 50, 1000) for Nikon 1026 plasmid\972_2000ms\Multichannel-0202.tif
mask_thres = 0.019;
area_low = 10;
area_high = 50;
area_hei = 1000;


%mask
%p = imread(['D:\Postdoc\microscope\0808\829+811\' filename]);

%turbo cells
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\turbo\962+935\Organoid2015.tif_Files\Organoid2015_c0.tif' filename]);
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\turbo\962+936 0 939ms\Organoid2010.tif_Files\Organoid2010_c0.tif' filename]);
%dh10b
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\Dh10b\962+936 0 550ms\Organoid1984.tif_Files\Organoid1984_c0.tif' filename]);
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\Dh10b\962+935\Organoid1987_110ms.tif_Files\Organoid1987_c0.tif' filename]);

% import data
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\nikon 1104\936 new 0\Multichannel-0702.tif' filename]);
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\Nikon 1026 plasmid\972_2000ms\Multichannel-0202.tif' filename]);
p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\nikon 1211 plasmid different copy number\871+963 induced 3hours 2000ms\Multichannel-test1.tif' filename]);
p2 = p;
pm = im2double(p);

%cell masks
m2 = zeros(size(p));
% identify all spots that have been collected
m3 = zeros(size(p));
m4 = zeros(size(p)*4);

%should be pm < thres for constrast and > thres for fluorescence
%original thres 0.02 for fluorescence images
%original thres 0.047 for phase contrast
%get mask from fluorescence images
%p_br2 = im2double(p_br);
m_br = double(pm>mask_thres);

%2017-12-12 should use a repetative background removing strategy
CC=bwconncomp(m_br);
stats=regionprops(CC,'basic');
num = 0;
for i = 1:CC.NumObjects
    larea(i) = stats(i).Area;
    if larea(i) > 100 && larea(i) < 4000000
        num = num + 1;
        m2(CC.PixelIdxList{i}) = num;
        p2(CC.PixelIdxList{i}) = p(CC.PixelIdxList{i})-prctile(p(CC.PixelIdxList{i}),60);
        %intensity(num) = sum(sum(p(CC.PixelIdxList{i})));
        %p2(CC.PixelIdxList{i}) = p(CC.PixelIdxList{i})-imgaussfilt(p(CC.PixelIdxList{i}),5,'FilterDomain','spatial');
        %cf(num,1) = (mean(p2(CC.PixelIdxList{i}))-mean(mean(p2(1:10,1:10))))/exp_t;
        %cf(num,2) = larea(i);
        
        %second round of background removal
        cella = zeros(size(p2));
        cella(CC.PixelIdxList{i}) = p2(CC.PixelIdxList{i});
        CC2=bwconncomp(cella>300);
        stats2=regionprops(CC2,'basic');
        for k = 1:CC2.NumObjects
            larea2(k) = stats2(k).Area;
            if larea2(k) > 200 && larea2(k) < 4000000
                cella(CC2.PixelIdxList{k}) = cella(CC2.PixelIdxList{k})-prctile(cella(CC2.PixelIdxList{k}),60);
            end
        end
        p2(CC.PixelIdxList{i}) = cella(CC.PixelIdxList{i});
    end
end


%subtracting background
%get spot intensity version 1
% m = double(pm>thres);
% %m = double(im2double(p2)>thres);
% CC=bwconncomp(m);
% stats=regionprops(CC,'basic');
% num = 0;
% %[s1,s2] = size(m4);
% for i = 1:CC.NumObjects
%     larea(i) = stats(i).Area;
%     %if larea(i) > 100 && larea(i) < 2000
%     %if larea(i)<50 && larea(i) >10
%     if larea(i)<area_high && larea(i) >area_low && max(max(p(CC.PixelIdxList{i})))> area_hei ...
%             && stats(i).Area/max(stats(i).BoundingBox(3:4))^2>0.5
%         num = num + 1;
%         m3(CC.PixelIdxList{i}) = num;
%         %m4(CC.PixelIdxList{i}*16) = num;
%         %yindex(num) = mod(CC.PixelIdxList{i}(1)*16-1,s1)+1;
%         %xindex(num) = floor(CC.PixelIdxList{i}(1)*16/s1);
%         %p2(CC.PixelIdxList{i}) = p(CC.PixelIdxList{i})-median(median(p(CC.PixelIdxList{i})));
%         %how to get the background signal for those pixels
%         intensity(num) = sum(sum(p2(CC.PixelIdxList{i})));
%         %intstr{num} = num2str(intensity(num));
%         %p2(CC.PixelIdxList{i}) = p(CC.PixelIdxList{i})-imgaussfilt(p(CC.PixelIdxList{i}),5,'FilterDomain','spatial');
%         %cf(num,1) = (mean(p2(CC.PixelIdxList{i}))-mean(mean(p2(1:10,1:10))))/exp_t;
%         %cf(num,2) = larea(i);
%     end
% end

%get spot intensity version 2
num = 0;
 for j = 1:max(max(m2))
        cella = zeros(size(p2));
        cella(find(m2==j)) = p2(find(m2==j));
        %parameter for 1105
        %cella = 2000 and max(max(cella)) = 4000 area (10-50)
        CC=bwconncomp(cella>300);
        stats=regionprops(CC,'basic');
        for k = 1:CC.NumObjects
            larea(k) = stats(k).Area;
            if larea(k)<500 && larea(k) >10 && max(max(cella))> 600 ...
                    && stats(k).Area/max(stats(k).BoundingBox(3:4))^2>0.5 ...
                    && 1
                num = num + 1;
                m3(CC.PixelIdxList{k}) = num;
                intensity(num) = sum(sum(cella(CC.PixelIdxList{k})));
            end
        end
 end

subplot(2,2,1)
image(m3*100)
 %hold on
 %text(xindex-5,yindex+10,intstr,'FontSize',0.2);
subplot(2,2,2)
%initial cell mask
image(m_br*100)
subplot(2,2,3)
%image(m*100)
subplot(2,2,4)
%fluorescence image after the removal of background
image(p2/10)
%hist(intensity)

filename = ['D:\Dropbox (MIT)\Postdoc\programs\Image Process\microscopy\' 'temp_proc.tif'];
imwrite(p2,filename);
%output files
% if nargin>2
%     load fishseg002.mat
%     LcFull = m2;
%     eval(['save fishseg' num2str(fishnum,'%03d') '.mat LN timestamp ph tempsegcorrect phaseFullSize Lc LNsub phsub rect LcFull;'])
%     copyfile(['fishseg' num2str(fishnum,'%03d') '.mat'],'D:\Postdoc\programs\spatzcells\Mysamples\segmentation\masks');
% 
%     %filename = ['D:\Postdoc\programs\spatzcells\Mysamples\images\fish_' num2str(sample,'%03d') 'xy' num2str(frame,'%02d') 'z1c2.tif'];
%     %imwrite(p2,filename);
% end


end

