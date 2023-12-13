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

[c1,l1] = wavedec2(Ared, 4, 'bior4.4');
[c2,l2] = wavedec2(Agreen, 4, 'bior4.4');
[c3,l3] = wavedec2(Ablue, 4, 'bior4.4');

% Compress by putting small coefficients to zero.
% Experiment with the threshold p.
T = max([max(abs(c1)), max(abs(c2)), max(abs(c3))]);
p = 1e-2;
threshold = p*T
I1 = find(abs(c1) < threshold);
c1(I1) = 0;
I2 = find(abs(c2) < threshold);
c2(I2) = 0;
I3 = find(abs(c3) < threshold);
c3(I3) = 0;

A2red = waverec2(c1, l1, 'bior4.4');
A2green = waverec2(c2, l1, 'bior4.4');
A2blue = waverec2(c3, l1, 'bior4.4');
A2 = zeros(size(A));
A2(:,:,1) = A2red;
A2(:,:,2) = A2green;
A2(:,:,3) = A2blue;
A2 = uint8(A2);

% Plot the result after compression
figure
imshow(A2)

% Compression ratio
N = length(c1)+length(c2)+length(c3);
M = length(I1)+length(I2)+length(I3);
compression_ratio = N / (N-M)

% Root mean square error
mserr = sqrt(sum(sum(sum( (A2-A).^2 ))) / length(A(:)))
% the norm of A(:) is equal to the so-called Frobenius norm
mserr_rel = mserr / norm(double(A), "fro")
