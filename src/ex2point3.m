% Read the image. In this case, this results in a
% M x N x 3 tensor, with RGB components.
% The colormap is empty.
[A,cmap] = imread('bib.png');
figure
imshow(A)

% The tensor contains red, green and blue components
Ared = A(:,:,1);
Agreen = A(:,:,2);
Ablue = A(:,:,3);

% Recent versions of the wavelet toolbox can transform the tensor A.
% It is generally safer to split the channels and do it manually
%[c,l] = wavedec2(A, 4, 'bior4.4');

[ca1,chd1,cvd1,cdd1] = swt2(Ared(:,1:1348), 2,'db6');
[ca2,chd2,cvd2,cdd2]= swt2(Agreen(:,1:1348), 2,'db6');
[ca3,chd3,cvd3,cdd3] = swt2(Ablue(:,1:1348),2,'db6');

