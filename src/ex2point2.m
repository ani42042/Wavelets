% Read the image. In this case, this results in a
% M x N x 3 tensor, with RGB components.
% The colormap is empty.
[A,cmap] = imread('bib.png');
figure
imshow(A)

%% compress image
p = 9e-3; thresholdType = 'Soft'; waveletType = 'db30';
A_noise = A;
[A2,mserr,mserr_rel,compression_ratio,threshold,SNR] = denoising_image(A_noise,p,thresholdType,waveletType,A);

% Plot the result after denosing/compression
figure
imshow(A2)

%% Add noise to image

% add noise to image
rng(42)

A_noise = imnoise(A,'salt & pepper',0.01);
%A_noise = imnoise(A,'speckle');
%A_noise = imnoise(A,"gaussian");
SNR_noise = 10*log10(norm(double(A),'fro')/(norm(double(A)-double(A_noise),'fro')));
msNoise = norm(double(A_noise)-double(A),'fro');
% Plot the noisy image
figure
imshow(A_noise)

%% Denoise image for specific parameters
p = 6e-3; thresholdType = 'Soft'; waveletType = 'db30'; 
[A2,mserr,mserr_rel,compression_ratio,threshold,SNR] = denoising_image(A_noise,p,thresholdType,waveletType,A);

% Plot the result after denosing/compression
figure
imshow(A2)

%% computing the best parameters for denosing this image
wavelets = ["bior1.1", "bior1.3", "bior1.5","bior2.2", "bior2.4",...
            "bior2.6", "bior2.8","bior3.1","bior3.3", "bior3.5",...
            "bior3.7","bior3.9", "bior4.4", "bior5.5", "bior6.8",...
            "rbio1.1", "rbio1.3", "rbio1.5","rbio2.2", "rbio2.4",...
            "rbio2.6", "rbio2.8","rbio3.1", "rbio3.3", "rbio3.5",...
            "rbio3.7","rbio3.9", "rbio4.4", "rbio5.5", "rbio6.8"];

tempStrings = strings(0);
for i = 1:45
    tempStrings(i) = strcat("db",num2str(i));
end

% the parameters to be tested
wavelets = [wavelets, tempStrings];
thesholding = ["Hard","Soft"];
p = linspace(1e-3,1,10);

% parameters to store results
SNRMat = zeros(length(thesholding),length(wavelets),length(p));
ThreshMat = zeros(length(thesholding),length(wavelets),length(p));
%%
n = length(wavelets)*length(thesholding); index = 1;
progressbar
for i = 1:length(thesholding)
    thesholdType = thesholding(i);
    for j = 1:length(wavelets)
        wave = wavelets(j);
        for k = 1:length(p)
            delta = p(k);
            [~,~,~,~,threshold,SNR] = denoising_image(A_noise,delta,thresholdType,wave,A);
            SNRMat(i,j,k) = SNR;
            ThreshMat(i,j,k) = threshold;
        end
        progressbar(index/n);
        index = index + 1;
    end
end
%save Data.mat ThreshMat SNRMat
%% plot best image
load Data.mat

ibest = 1;
jbest = 1;
kbest = 1;
bestSNR = SNRMat(ibest,jbest,kbest);
for i = 1:size(SNRMat,1)
    for j=1:size(SNRMat,2)
        for k=1:size(SNRMat,3)
            if (bestSNR <= SNRMat(i,j,k))
                bestSNR = SNRMat(i,j,k); ibest = i; jbest = j; kbest = k;
            end
        end
    end
end
bestSNR
BestThresh = ThreshMat(ibest,jbest,kbest);
[BestImage,~,~,~,TR,SNR] = denoising_image(A_noise,p(kbest),thesholding(ibest),wavelets(jbest),A);
% Plot the result after denosing for best params
figure
imshow(BestImage)

% function to apply hard thresholding
function [c,I] = Hard_threshold(delta, c)
    I = find(abs(c) < delta);
    c(I) = 0;
end

% function to apply soft thresholding
function [c, I] = Soft_threshold(delta,c)
    I = find(abs(c) < delta);
    c = sign(c).*(abs(c)-delta);
    c(I) = 0;
end

% function to desnoise an image
function [A2,mserr,mserr_rel,compression_ratio,threshold,SNR] = denoising_image(A,p,threshtype,waveletType,A_correct)
    % if there is no A_correct, interprete as compression and not
    % denosing
    if nargin < 5
        A_correct = A;
    end

    % The tensor contains red, green and blue components
    Ared = A(:,:,1);
    Agreen = A(:,:,2);
    Ablue = A(:,:,3);
    
    % Recent versions of the wavelet toolbox can transform the tensor A.
    % It is generally safer to split the channels and do it manually
    %[c,l] = wavedec2(A, 4, 'bior4.4');
    
    [c1,l1] = wavedec2(Ared, 4, waveletType);
    [c2,l2] = wavedec2(Agreen, 4, waveletType);
    [c3,l3] = wavedec2(Ablue, 4, waveletType);
    
    % Compress by putting small coefficients to zero.
    % Experiment with the threshold p.
    T = max([max(abs(c1)), max(abs(c2)), max(abs(c3))]);
    %p = 1e-1;
    threshold = p*T;
    if (strcmp("Hard",threshtype))
        [c1,I1] = Hard_threshold(threshold, c1);
        [c2,I2] = Hard_threshold(threshold, c2);
        [c3,I3] = Hard_threshold(threshold, c3);
    else
        [c1,I1] = Soft_threshold(threshold, c1);
        [c2,I2] = Soft_threshold(threshold, c2);
        [c3,I3] = Soft_threshold(threshold, c3);
    end
    A2red = waverec2(c1, l1, waveletType);
    A2green = waverec2(c2, l2, waveletType);
    A2blue = waverec2(c3, l3, waveletType);
    A2 = zeros(size(A));
    A2(:,:,1) = A2red;
    A2(:,:,2) = A2green;
    A2(:,:,3) = A2blue;

    % if condition true we are denoising, thus can compute the signal to
    % noise ratio
    if nargin >= 5
        SNR = 10*log10(norm(double(A_correct),'fro')/(norm(double(A_correct)-A2,'fro')));
    end
    
    A2 = uint8(A2);
    
    % Plot the result after compression
    % figure
    % imshow(A2)
    
    % Compression ratio
    N = length(c1)+length(c2)+length(c3);
    M = length(I1)+length(I2)+length(I3);
    compression_ratio = N / (N-M);
    
    % error
    mserr = norm(double(A2)-double(A_correct),'fro');
    mserr_rel = mserr / norm(double(A_correct), "fro");
end