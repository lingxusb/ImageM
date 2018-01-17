function [ ppad ] = padding( p2 )
%PADDING padding the existing image p2 by box_d
%   the same size of image is copied to padding region
%   2018-1-12
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
end

