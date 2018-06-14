%FLURTOPHASE convert fluorescence image to phase images for schneitz cells
%   2018-01-11

%% specify all the paths and initialize schnitzcells
exp_date = '2018-06-08'
schn_path = 'D:\Dropbox (MIT)\Postdoc\programs\Schnitzcells\samples\';
p = initschnitz('TestSchnitz-01',exp_date,'e.coli',...
'rootDir',schn_path);
source_dir = 'D:\Dropbox (MIT)\Postdoc\microscope\nikon 20180608 copy number 100X\877 g 200ms r 2s\';
save_dir = [schn_path exp_date '\TestSchnitz-01\images\'];

filename{1} = 'Multichannel-0003.tif';
filename{2} = 'Multichannel-0203.tif';
filename{3} = 'Multichannel-0403.tif';
filename{4} = 'Multichannel-0403.tif';
filename{5} = 'Multichannel-0403.tif';
filename{6} = 'Multichannel-0403.tif';

crop = 1;

if crop == 1
    crop_image = imread('D:\Dropbox (MIT)\Postdoc\programs\Image Process\GitRep\crop_0604.tif');
end


%% process all the images for cell segmentation

for i = 1:6
    imgr = imread([source_dir filename{i}]);

    
    %get the region for cells
    mthres = [1100,100];
    pm = double(imgr);
    spots = zeros(size(pm));
    
    % Remove background and get spot intensity using the old methods.
    %[pm,spots] = delBack(pm, mthres,0, zeros(size(pm)));
    
    % using imfindcircles to find spots
    % current parameters: radii from 1 to 4, sensitivity 0.98-1, should get
    % rid of the bright cells at the first time or drop some spots in bright
    % cells afterwards.
    % TwoStage seems to be better than the default method. Sensitivity can
    % be set to 1
    %image(pm/100)
    %[centersBright, radiiBright] = imfindcircles(pm,[1 4],'ObjectPolarity','bright','Sensitivity',1,'Method','TwoStage');
    %viscircles(centersBright, radiiBright,'Color','b');
    
%   remove over expression ONLY for some of the images
%     filen = filename{i};
%     imgr2 = imread([source_dir [filen(1:length(filen)-6) num2str(1,'%02d') '.tif']]);
%     m_br = double(imgr2==16383);
%     CC=bwconncomp(m_br);
%     stats=regionprops(CC,'basic');
%     for j = 1:1:CC.NumObjects
%         larea(j) = stats(j).Area;
%         if larea(j)>150
%             imgr(CC.PixelIdxList{j}) = 740;
%         end
%     end

    %crop image only for some of the images
    if crop == 1
        imgr(find(crop_image<max(max(crop_image)))) = 600;
    end
    
%   Fluore mask only
    % prefix -01-t- indicates that we used these files for segmentation
    % based on fluorescence images
    write_name = [schn_path...
      exp_date '\TestSchnitz-01\images\TestSchnitz-01-t-' num2str(i,'%03d') '.tif'];
    imwrite(imgr,write_name);

%   PHASE constrast mask only
    % based on phase images
%         write_name = [schn_path...
%          exp_date '\TestSchnitz-01\images\TestSchnitz-01-p-' num2str(i,'%03d') '.tif'];
%     imwrite(imgr,write_name);

    %write_name = ['D:\Dropbox (MIT)\Postdoc\programs\Schnitzcells\samples\'...
    %     exp_date '\TestSchnitz-01\images\TestSchnitz-01-spots-t-' num2str(i,'%03d') '.tif'];
    %imwrite(uint16(spots),write_name);
end

%   Fluore mask only
p = segmoviefluor(p,'maxThresh',0.1);

%   PHASE constrast mask only
%p = segmoviephase(p,'maxThresh',0.1,'minThresh',0.1);

p = manualcheckseg(p);

%load masks
for i = 1:6
    seg_path = load([schn_path exp_date '\TestSchnitz-01\segmentation\TestSchnitz-01seg' num2str(i,'%03d') '.mat'],'Lc');
    mask{i} = seg_path.Lc;
end

%% remove background from samples and identify spots (old method)
%   pm: gray-scale image
%   mthres: threshold to identify cells, thres1 in the first round and thres2 in
%   the following rounds 
%   iter: record the iteration depth
%   spots: spots map
%   2018-01-12
function [pm,spots] = delBack(pm, mthres, iter,spots)
    %mask and fill all the regions
    if iter == 0
        m_br = double(pm>mthres(1));
    else
        m_br = double(pm>mthres(2));
    end
    CC=bwconncomp(m_br);
    stats=regionprops(CC,'basic');
    for j = 1:CC.NumObjects
        larea(j) = stats(j).Area;
        %identify all the cells and remove background
        if larea(j) > 100 && larea(j) < 4000000
            pm(CC.PixelIdxList{j}) = pm(CC.PixelIdxList{j})-prctile(pm(CC.PixelIdxList{j}),60);
            % function recursion
            if iter < 2
                cella = zeros(size(pm));
                cella(CC.PixelIdxList{j}) = pm(CC.PixelIdxList{j});
                [cella, spots] = delBack(cella, mthres, iter+1,spots);
                pm(CC.PixelIdxList{j}) = cella(CC.PixelIdxList{j});
            end        
        end
        %identfiy spots
        if iter > 0
            if larea(j)<500 && larea(j) >10 && max(max(pm(CC.PixelIdxList{j})))> 300 ...
                    && stats(j).Area/max(stats(j).BoundingBox(3:4))^2>0.5
                spots(CC.PixelIdxList{j}) = pm(CC.PixelIdxList{j});
            end
        end
    end
end

%% pre-process image 
%   should make bright cells darker and ehance constrast.
function pm = delBrig(pm, thres)
%     m_br = double(pm>thres);
%     CC=bwconncomp(m_br);
%     stats=regionprops(CC,'basic');
%     for j = 1:CC.NumObjects
%         larea(j) = stats(j).Area;
%         %identify all the cells and remove background
%         if larea(j) > 100 && larea(j) < 4000000
%             pm(CC.PixelIdxList{j}) = pm(CC.PixelIdxList{j})-prctile(pm(CC.PixelIdxList{j}),50)+2100;
%         end
%     end
end

