function [ spots,rfp, gfp, background, intensity, mintensity ] = MasktoSpots( mask,filename, schn_path, exp_date, index)
%MASKTOSPOTS generate the statistic of spots in a mask
%   2018-01-17
%   spots: the intensity of spots
%   rfp: total spots intensity per cell
%   gfp: fluorescence for background
%   background: background for spot channel (substracted)
%   intensity: total fluorescence of spot channel for each cell
%   mintensity: mintensity for each cell

%read the fluorescence file 
imgr = imread([filename(1:length(filename)-6) num2str(1,'%02d') '.tif']);
s = size(imgr);

%   PHASE constrast mask only
%   crop fluorescence image when use phase constrast to mask (size doesn't fit)
% imgr(s(1)-5:s(1),:) = [];
% imgr(1:6,:) = [];
% imgr(:,s(2)-3:s(2)) = [];
% imgr(:,1:5) = [];
% s = size(imgr);

% import gfp channel
imgr2 = imread([filename]);

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
background = zeros(max2(mask{index}),1);
for i = 1:max2(mask{index})
    maski = delout(find(mask{index}==i)+drift, s(1)*s(2));
    background(i) = prctile(imgr(maski),60);
    background_g(i) = double(prctile(imgr2(maski-drift),60));
    intensity(i) = sum(sum(imgr(maski)));
    mintensity(i) = sum(sum(imgr(maski)))/sum(sum(mask{index}==i));
    %   Fluor mask only
    background_all = uint16(background(i)*double(imgr2(maski-drift))/background_g(i));
    background_all(background_all< prctile(imgr(maski),60)) = prctile(imgr(maski),60);
    imgr(maski) = minus(imgr(maski),background_all);
%     cella = zeros(size(imgr));
%     cella(maski) =  imgr(maski);
%     PSF = fspecial('gaussian',5,1);
%     cella = imfilter(cella,PSF,'symmetric','conv');
%     m_br = double(cella<500).*double(cella>0);
%     CC=bwconncomp(m_br);
%     stats=regionprops(CC,'basic');
%     for j = 1:1:CC.NumObjects
%         larea(j) = stats(j).Area;
%         if stats(j).Area/max(stats(j).BoundingBox(3:4))^2 < 0.4
%             imgr(CC.PixelIdxList{j}) = 0;
%         end
%     end
    %   PHASE constrast mask only
    %imgr(maski) = minus(imgr(maski),background(i));
end


%% expand mask
for i = 1:max2(mask{index})
    maski = delout(find(mask{index}==i)+drift, s(1)*s(2));
    for j = 1:0
        additional = [maski-s(1),maski+s(1),maski-1,maski+1];
        additional = reshape(additional,4*length(maski),1);
        additional(additional<=0) = [];
        additional(additional>=s(1)*s(2)) = [];
        additional = unique(additional);
        additional_pos = additional-drift;
        additional_pos(additional_pos<=0) = [];
        additional_pos(additional_pos>s(1)*s(2)) = [];
        additional(find(mask{index}(additional_pos)>0)) = [];
        additional = setdiff(additional,maski);
        if sum(imgr(additional)>background(i)) <2
            break;
        end
        %choosing this condition depending on what you want to see, a fully
        %expansion or only the pixels that make sense.
        add_element = ones(size(additional));%(imgr(additional)>background(i));%.*(imgr(additional)< max(imgr(maski)));
        add_element = find(add_element == 1);
        maski = [maski; additional(add_element)];
        imgr(additional(add_element)) = minus(imgr(additional(add_element)),background(i));
    end
    mask_after{i} = maski;
end

%% calculate statistics

% get gfp statistics
for i = 1:max2(mask{index})
    gfp(i) = prctile(imgr2(find(mask{index}==i)),50);
end

%get spots intensity
rfp =  zeros(max2(mask{index}),1);
img_spots = zeros(size(imgr));
img_spots_mask = zeros(size(imgr));
spots = [];
s = size(imgr);
for i = 1:max2(mask{index})
    cella = zeros(size(imgr));
    %maski = delout(find(mask{index}==i)+drift, s(1)*s(2));
    maski = mask_after{i};
    cella(maski) = imgr(maski);
    %filter the fluorescence image
    %original fspecial('gaussian',10,1)
    %cella(cella < 600) = 0;
    PSF = fspecial('gaussian',5,1);
    Filtered = imfilter(cella,PSF,'symmetric','conv');
    %Filtered = cella;
    imgr(maski) =  Filtered(maski);
    box_size = 11;
    [i_map, j_map] = meshgrid(1:s(2),1:s(1));
    all_maxima = imregionalmax(Filtered, 8);
    maxima_i = find(all_maxima==1);
    maxima_i = remove_spot_boundary(maxima_i, maski, s);
    for j = 1:length(maxima_i)
        pos(1)= mod(maxima_i(j)-1,s(1)) + 1;
        pos(2)= floor(maxima_i(j)/s(1));
        if Filtered(maxima_i(j))> 400 && pos(1) > floor(box_size/2) && pos(2) > floor(box_size/2) &&...
                pos(1) < s(1) - floor(box_size/2) && pos(2) < s(2) - floor(box_size/2)
            Pixels_to_analyze = zeros(size(imgr));
            Pixels_to_analyze(pos(1)-floor(box_size/2):pos(1)+floor(box_size/2),pos(2)-floor(box_size/2):pos(2)+floor(box_size/2)) = 1;
            Pixels_to_analyze= find(Pixels_to_analyze==1);
            Pixels_to_analyze = intersect(Pixels_to_analyze,maski);
            if max2(img_spots_mask(Pixels_to_analyze)==i) == 0
                X1 = i_map(Pixels_to_analyze);
                X2 = j_map(Pixels_to_analyze);
                Y = imgr(Pixels_to_analyze);
                spot_height = double(Filtered(maxima_i(j)));
                p0 = [spot_height,pos(1),pos(2),1,1,0,0];
                lb = [0,pos(1)-5,pos(2)-5,0.01* 10000 / spot_height,0.01* 10000 / spot_height,0,0];
                ub = [double(Filtered(maxima_i(j)))+1000,pos(1)+5,pos(2)+5,10,10,0.1,1000];
                [p, Y_pred, resnorm] = fitGuassian(X1, X2, Y, p0, lb, ub, pos);
                %maximum in the cell mask -- may not need this condition
                %if sum(maski == (floor(p(3))-1)*s(1)+floor(p(2))) > 0
                spots = [spots;cat(2,p,resnorm)];
                rfp(i) = rfp(i) + p(1);
                img_spots(Pixels_to_analyze) = Y_pred;
                img_spots_mask(Pixels_to_analyze) = i;
                %end
            end
        end
    end
    
%     m_br = double(cella>600);
%     CC=bwconncomp(m_br);
%     stats=regionprops(CC,'basic');
%     for j = 1:1:CC.NumObjects
%         larea(j) = stats(j).Area;
%         if larea(j)<100 && larea(j) >10 && max(max(cella(CC.PixelIdxList{j})))> 1200 ...
%                     && stats(j).Area/max(stats(j).BoundingBox(3:4))^2>0.3
%                 spots = [spots, sum(sum(cella(CC.PixelIdxList{j})))];
%                 rfp(i) = rfp(i) + sum(sum(cella(CC.PixelIdxList{j})));
%                 if sum(sum(cella(CC.PixelIdxList{j}))) < 20000
%                     img_spots(CC.PixelIdxList{j}) = cella(CC.PixelIdxList{j});
%                 end
%         end
%     end
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

%% fit guassian distribution
function [p, Y_pred, resnorm] = fitGuassian(X1, X2, Y, p0, lb, ub, pos)
gau2_gen = @(n) ['p(' num2str(1+6*(n-1)) ').*', ... %%% 2D single-gaussian model function text generator
                 'exp(', ...                      %%% exp(-( a(x-x0)^2 + b(y-y0)^2 + 2c(x-x0)(y-y0) ))
                 '-(', ...
                 '(p(' num2str(4+6*(n-1)) ')).*(x(:,1)-p(' num2str(2+6*(n-1)) ')).^2+', ...
                 '(p(' num2str(5+6*(n-1)) ')).*(x(:,2)-p(' num2str(3+6*(n-1)) ')).^2+', ...
               '2*(p(' num2str(6+6*(n-1)) ')).*(x(:,1)-p(' num2str(2+6*(n-1)) ')).*(x(:,2)-p(' num2str(3+6*(n-1)) '))', ...
                 '))','+p(' num2str(7+6*(n-1)) ')'];  
eval(['f = @(p,x)' gau2_gen(1) ';']);
options = optimset('Display','off');
[p, resnorm, ~, exitflag] = lsqcurvefit(f, p0, double([X2 X1]), double(Y), lb, ub,options);
a = p(4);
b = p(5);
c = p(6);
if a*b - c^2 < 0 || c > 0.1% || min(a,b) < 0.01* 10000 / p(1) %a*b-c^2 < 0 || abs(a-b)>1 || c > 0.3 || min(a,b) < 0.01* 20000 / p(1)
    'spot is poorly fitted'
    Y_pred = 0;
    p = zeros(1,8);
else
    spot_intensity = p(1)*pi/sqrt(a.*b-c.^2);
    Y_pred = f(p,[X2 X1])- p(7);
    p = [spot_intensity p];
end
end


%% remove maxima near the boundary
function maxima_i = remove_spot_boundary(maxima_i, maski, s)
     for i = 1:length(maxima_i)
%          t(i) = length(intersect([maxima_i(i) + 1, maxima_i(i) - 1, maxima_i(i) + s(1), maxima_i(i) - s(1),...
%              maxima_i(i) + 1 + s(1), maxima_i(i) - 1 + s(1), maxima_i(i) + 1 - s(1), maxima_i(i) - 1 - s(1)],maski));
         t(i) = length(intersect([maxima_i(i) + 1, maxima_i(i) - 1, maxima_i(i) + s(1), maxima_i(i) - s(1)],maski));
     end
     maxima_i(find(t<4)) = [];
end
