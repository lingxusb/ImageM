function [ spots, gfp ] = MasktoSpots( mask,filename, schn_path, exp_date, index)
%MASKTOSPOTS generate the statistic of spots in a mask
%   2018-01-17

%read the fluorescence file 
imgr = imread([filename(1:length(filename)-6) num2str(2,'%02d') '.tif']);
s = size(imgr);

%% calculate drif of images
mask2 = mask{index};
mask2(1:8,:) = 0;
mask2(:,1:8) = 0;
mask2(s(1)-8:s(1),:) = 0;
mask2(:,s(2)-8:s(2)) = 0;
for i = 1:16
    for j = 1:16
        inten(i,j) = sum(sum(imgr(find(mask2>0)+(i-8)*s(1)+j-8)));
    end
end
minten = max(inten(:));
[rindex, cindex] = find(inten == minten);
drift = (rindex-8)*s(1)+cindex-8;

%% remove background
for i = 1:max2(mask{index})
    maski = delout(find(mask{index}==i)+drift, s(1)*s(2));
    background = prctile(imgr(maski),60);
    imgr(maski) = imgr(maski)-background;
end

%% calculate statistics
% import gfp channel
imgr2 = imread([filename]);
% get gfp statistics
for i = 1:max2(mask{index})
    gfp(i) = prctile(imgr2(find(mask{index}==i)),50);
end

%get spots intensity
spots =  zeros(max2(mask{index}),1);
img_spots = zeros(size(imgr));
for i = 1:max2(mask{index})
    cella = zeros(size(imgr));
    maski = delout(find(mask{index}==i)+drift, s(1)*s(2));
    cella(maski) = imgr(maski);
    m_br = double(cella>200);
    CC=bwconncomp(m_br);
    stats=regionprops(CC,'basic');
    for j = 1:1:CC.NumObjects
        larea(j) = stats(j).Area;
        if larea(j)<150 && larea(j) >10 && max(max(cella(CC.PixelIdxList{j})))> 300 ...
                    && stats(j).Area/max(stats(j).BoundingBox(3:4))^2>0.4
                spots(i) = spots(i) + sum(sum(cella(CC.PixelIdxList{j})));
                img_spots(CC.PixelIdxList{j}) = cella(CC.PixelIdxList{j});
        end
    end
end

%% output the image that removed the background
eval(['mkdir ''',[schn_path exp_date '\TestSchnitz-01\'],''' ','spots'])
write_name = [schn_path...
    exp_date '\TestSchnitz-01\spots\TestSchnitz-01-t-' num2str(index,'%03d') '.tif'];
imwrite(imgr,write_name);

write_name = [schn_path...
    exp_date '\TestSchnitz-01\spots\TestSchnitz-01-t-spots' num2str(index,'%03d') '.tif'];
imwrite(uint16(img_spots),write_name);

end

%% remove outlier of the drift
% input arg: index is for mask of cell i, s is the size of the image.
function index = delout(index, s)
index(index<0) = [];
index(index > s) = [];
end

