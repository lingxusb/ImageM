%FLURTOPHASE convert fluorescence image to phase images for schneitz cells
%   2018-01-11

%% specify all the paths and initialize schnitzcells
exp_date = '2018-01-11'
p = initschnitz('TestSchnitz-01',exp_date,'e.coli',...
'rootDir','D:\Dropbox (MIT)\Postdoc\programs\Schnitzcells\samples');
source_dir = 'D:\Dropbox (MIT)\Postdoc\microscope\nikon 1211 plasmid different copy number\871+963 induced 3hours 2000ms\';
save_dir = ['D:\Dropbox (MIT)\Postdoc\programs\Schnitzcells\samples\' exp_date '\TestSchnitz-01\images\'];

filename{1} = 'MultichannelF-t1.tif';
filename{2} = 'MultichannelF-t2.tif';
filename{3} = 'MultichannelF3.tif';
filename{4} = 'MultichannelF4.tif';
filename{5} = 'MultichannelF5.tif';
filename{6} = 'MultichannelF.tif';

%% process all the images for cell segmentation
for i = 1:1
    imgr = imread([source_dir filename{i}]);
    p2 = imgr;
    %2100 is a good threshold
    %p2(p2>2100)  = 2100;

    %prefix -01-t- indicates that we used these files for segmentation
    %based on fluorescence images
    
    %scaffold
    box_d = 25; 
    img_d = size(p2);
    ppad = zeros(size(p2) + 2*box_d,'uint16');
    ppad(box_d+1:box_d+img_d(1),box_d+1:box_d+img_d(2)) = p2;
    %padding
    for j = 1:img_d(2)-box_d+1
        ppad(1:box_d, j+box_d:j+2*box_d-1) = ppad(box_d+1:box_d*2, j+box_d:j+2*box_d-1);
        ppad(img_d(1)+1+box_d:img_d(1)+2*box_d, j+box_d:j+2*box_d-1) = ppad(img_d(1)+1:img_d(1)+box_d, j+box_d:j+2*box_d-1);
    end
    for j = 1:img_d(1)+box_d+1
        ppad(j:j+box_d-1,1:box_d) = ppad(j:j+box_d-1, box_d+1:box_d*2);
        ppad(j:j+box_d-1,img_d(2)+box_d+1:img_d(2)+2*box_d) = ppad(j:j+box_d-1, img_d(2)+1:img_d(2)+1*box_d);
    end
    
    for m = 1:img_d(1)
        for n = 1:img_d(2)
            window = ppad(m+1+box_d:m+2*box_d, n+1+box_d:n+2*box_d);
            p2(m,n) = p2(m,n) - prctile(double(window(:)),60);
        end
    end
    x = box_d;
    y = box_d*6;
    window = ppad(x+1:x+box_d, y+1:y+box_d);
    window = window(:);
    %hist(log(double(window)),100)
    %prctile(log(double(window)),60)
        
    
    write_name = ['D:\Dropbox (MIT)\Postdoc\programs\Schnitzcells\samples\'...
         exp_date '\TestSchnitz-01\images\TestSchnitz-01-t-' num2str(i,'%03d') '.tif'];
    imwrite(p2,write_name);
end

%p = segmoviefluor(p);


