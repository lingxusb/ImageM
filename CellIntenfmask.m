function [intensity, cnum] = CellIntenfmask( filename)
%IMAGEPROC get intensity distribution for cells
%   2017-11-03

%   image processing is divided into two steps, in the first step,
%   background fluorescence is subtracted from each cell. In the second
%   step, the spot instensity is collected

%mask
%load('D:\Dropbox (MIT)\Postdoc\programs\Schnitzcells\samples\2017-11-12\TestSchnitz-01\segmentation\TestSchnitz-01seg001.mat');
%LNsub = AmplifyMask(LNsub,3);
%turbo cells
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\turbo\962+935\Organoid2015.tif_Files\Organoid2015_c0.tif' filename]);
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\turbo\962+936 0 939ms\Organoid2010.tif_Files\Organoid2010_c0.tif' filename]);
%dh10b
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\Dh10b\962+936 0 550ms\Organoid1984.tif_Files\Organoid1984_c0.tif' filename]);
%p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\0928 strains\Dh10b\962+935\Organoid1987_110ms.tif_Files\Organoid1987_c0.tif' filename]);
p = imread(['D:\Dropbox (MIT)\Postdoc\microscope\nikon 1104\935 new\Multichannel-1104_05.tif' filename]);

%mask p
CC=bwconncomp(p>1100);
stats=regionprops(CC,'basic');
m2 = zeros(size(p));
p2 = p;
num = 0;

max_area = 1500;
min_area = 100;
[s1,s2] = size(p2);
for i = 1:CC.NumObjects
    larea(i) = stats(i).Area;
    if larea(i) > min_area && larea(i) < max_area
        num = num + 1;
        m2(CC.PixelIdxList{i}) = num;
        p2(CC.PixelIdxList{i}) = p(CC.PixelIdxList{i});%-median(median(p(CC.PixelIdxList{i})));
        yindex(num) = mod(CC.PixelIdxList{i}(1)-1,s1)+1;
        xindex(num) = floor(CC.PixelIdxList{i}(1)/s1);
        %intensity(num) = sum(sum(p(CC.PixelIdxList{i})));
        %p2(CC.PixelIdxList{i}) = p(CC.PixelIdxList{i})-imgaussfilt(p(CC.PixelIdxList{i}),5,'FilterDomain','spatial');
        %cf(num,1) = (mean(p2(CC.PixelIdxList{i}))-mean(mean(p2(1:10,1:10))))/exp_t;
        %cf(num,2) = larea(i);
    end
end

p2 = p;
p3 = p;
cell_num = 0;
num = 0;

m3 = zeros(size(p));
p3 = p3*0;
LNsub = m2;
for i = 1:max(max(LNsub))
    if sum(sum(LNsub==i))>min_area && sum(sum(LNsub==i))<max_area
        num = num+1;
        thres = prctile(p2(find(LNsub==i)),25);
        p_temp = p2(find(LNsub==i))-thres;
        if sum(p_temp>200) > sum(LNsub==i)*0.5
            thres = prctile(p2(find(LNsub==i)),35);
            %thres = thres +1000;
        end
        %if thres >3000
        %    thres = 3000;
        %end
        p3(find(LNsub==i)) = p2(find(LNsub==i))-thres;
        intensity(num) = (sum(sum(p2(find(LNsub==i))-thres)))*6;
        cnum(num) = i;
        intstr{num} = num2str(intensity(num));
    end
end
%image(m3);
imwrite(p3,'temp_10.tif');
image(LNsub)
hold on
text(xindex-5,yindex+10,intstr,'FontSize',8);
imwrite(LNsub,'temp_mask_10.tif');
%subtracting background
% area_low = 10;
% area_high = 200;
% area_hei = 200;
% m = double(im2double(p)>thres);
% m3 = zeros(size(m));
% %m = double(im2double(p2)>thres);
% CC=bwconncomp(m);
% stats=regionprops(CC,'basic');
% num = 0;
% for i = 1:CC.NumObjects
%     larea(i) = stats(i).Area;
%     ind = LNsub(CC.PixelIdxList{i});
%     ind = ind(ind~=0);
%     cell_index = mode(ind);
%     if cell_index >0
%         if larea(i)<area_high && larea(i) >area_low && max(max(p(CC.PixelIdxList{i})))> area_hei
%             num = num + 1;
%             m3(CC.PixelIdxList{i}) = num;
%             intensity(num) = (sum(sum(p2(CC.PixelIdxList{i})-prctile(p2(find(LNsub==cell_index)),25))))*5;
%         end
%     end
% end
% image(m3);

% this strategy would seperate the spots with different masks.
% p3 = zeros(size(p));
% for j = 1:max(max(LNsub))
%         p2(find(LNsub==j)) = p(find(LNsub==j)) - prctile(p(find(LNsub==j)),25);
%         cella = zeros(size(p2));
%         cella(find(LNsub==j)) = p2(find(LNsub==j));
%         CC=bwconncomp(cella>500);
%         stats=regionprops(CC,'basic');
%         c_inten = 0;
%         for k = 1:CC.NumObjects
%             larea(k) = stats(k).Area;
%             if larea(k)<area_high && larea(k) >area_low && max(max(cella))> area_hei
%                 num = num + 1;
%                 p3(CC.PixelIdxList{k}) = num;
%                 intensity(num) = sum(sum(cella(CC.PixelIdxList{k})))*5;
%                 c_inten =  c_inten + sum(sum(cella(CC.PixelIdxList{k})))*5;
%             end
%         end
%         if sum(sum(LNsub==j))>100 && 1
%             cell_num = cell_num+1;
%             c_dis(cell_num) = c_inten;
%         end 
% end
% image(p3)
end





