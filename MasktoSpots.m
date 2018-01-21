function [ spots, gfp ] = MasktoSpots( mask,filename, schn_path, exp_date)
%MASKTOSPOTS generate the statistic of spots in a mask
%   2018-01-17

%read the fluorescence file 
imgr = imread([filename num2str(2,'%02d') '.tif']);
s = size(imgr);

%% calculate drif of images
for i = 1:16
    for j = 1:16
        inten(i,j) = sum(sum(imgr(find(mask>0)+(i-8)*s(1)+j-8)));
    end
end
minten = max(inten(:));
[rindex, cindex] = find(inten == minten);
drift = (rindex-8)*s(1)+cindex-8;

%% remove background
for i = 1:max2(mask)
    background = prctile(imgr(find(mask==i)+drift),50);
    imgr(find(mask==i)+drift) = imgr(find(mask==i)+drift)-background;
end

%% calculate statistics
% import gfp channel
imgr2 = imread([filename num2str(3,'%02d') '.tif']);
% get gfp statistics
for i = 1:max2(mask)
    gfp(i) = prctile(imgr2(find(mask==i)),50);
end

%get spots intensity
spots =  zeros(max2(mask),1);
img_spots = zeros(size(imgr));
for i = 1:max2(mask)
    cella = zeros(size(imgr));
    cella(find(mask==i)+drift) = imgr(find(mask==i)+drift);
    m_br = double(cella>200);
    CC=bwconncomp(m_br);
    stats=regionprops(CC,'basic');
    for j = 1:1:CC.NumObjects
        larea(j) = stats(j).Area;
        if larea(j)<150 && larea(j) >10 && max(max(cella(CC.PixelIdxList{j})))> 300 ...
                    && 1 %stats(j).Area/max(stats(j).BoundingBox(3:4))^2>0.5
                spots(i) = spots(i) + sum(sum(cella(CC.PixelIdxList{j})));
                img_spots(CC.PixelIdxList{j}) = cella(CC.PixelIdxList{j});
        end
    end
end

%% output the image that removed the background
eval(['mkdir ''',[schn_path exp_date '\TestSchnitz-01\'],''' ','spots'])
write_name = [schn_path...
    exp_date '\TestSchnitz-01\spots\TestSchnitz-01-t-' num2str(1,'%03d') '.tif'];
imwrite(imgr,write_name);

write_name = [schn_path...
    exp_date '\TestSchnitz-01\spots\TestSchnitz-01-t-spots' num2str(1,'%03d') '.tif'];
imwrite(uint16(img_spots),write_name);

end

