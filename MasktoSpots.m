function [  ] = MasktoSpots( mask,filename, schn_path, exp_date)
%MASKTOSPOTS generate the statistic of spots in a mask
%   2018-01-17

%read the fluorescence file
imgr = imread(filename);
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

%% output the image that removed the background
eval(['mkdir ''',[schn_path exp_date '\TestSchnitz-01\'],''' ','spots'])
write_name = [schn_path...
    exp_date '\TestSchnitz-01\spots\TestSchnitz-01-t-' num2str(1,'%03d') '.tif'];
imwrite(imgr,write_name);

end

