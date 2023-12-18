% Read the image. In this case, this results in a
% M x N x 3 tensor, with RGB components.
% The colormap is empty.
[A,cmap] = imread('bib.png');
figure
imshow(A)

%% Add noise to image

% add noise to image
rng(42)

%A_noise = imnoise(A,'salt & pepper',0.01);
%A_noise = imnoise(A,'speckle');
A_noise = imnoise(A,"gaussian");
SNR = 10*log10(norm(double(A_noise),'fro')/(norm(double(A_noise)-double(A),'fro')));
%msNoise = sqrt(sum(sum(sum( (A_noise-A).^2 ))) / length(A(:)));
msNoise = norm(double(A_noise)-double(A),'fro');
% Plot the noisy image
figure
imshow(A_noise)

%%
% Recent versions of the wavelet toolbox can transform the tensor A.
% It is generally safer to split the channels and do it manually
%[c,l] = wavedec2(A, 4, 'bior4.4');

%[ca1,chd1,cvd1,cdd1] = swt2(Ared(:,1:1348), 2,'db6');
%[ca2,chd2,cvd2,cdd2]= swt2(Agreen(:,1:1348), 2,'db6');
%[ca3,chd3,cvd3,cdd3] = swt2(Ablue(:,1:1348),2,'db6');

% make sure the image is rows and colums can be divided by 2^Level
Level = 2;
select = 1:1348; % hardcoded

A = A(:,select,:);
A_noise = A_noise(:,select,:);
Ared = A_noise(:,:,1);
Agreen = A_noise(:,:,2);
Ablue = A_noise(:,:,3);

figure
imshow(A_noise)

%%

% parameters for decomposition
wave = "db1"; rel_threhold = 1e-2;
type = "redundant";
%type = "other";
threshold_type = "Hard";

[A2red,A2green,A2blue,temp] = denoise_tensor(Ared,Agreen,Ablue,wave,rel_threhold,type,threshold_type,Level);

A2 = zeros(size(A));
A2(:,:,1) = A2red;
A2(:,:,2) = A2green;
A2(:,:,3) = A2blue;
A2 = uint8(A2);

SNR_denoise = 10*log10(norm(double(A2),'fro')/(norm(double(A2)-double(A),'fro')));

% Plot the result after compression
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
Level = 2;
type = "standard";
TypeTransform = ["redundant","standard"];
% parameters to store results
SNRMat = zeros(length(TypeTransform),length(thesholding),length(wavelets),length(p));
ThreshMat = zeros(length(TypeTransform),length(thesholding),length(wavelets),length(p));
% mprev = +inf;
% BestImage = A_noise;
progressbar
index = 1; n = length(wavelets)*length(thesholding)*length(TypeTransform);
for l = 1:length(TypeTransform)
    type = TypeTransform(l);
    for i = 1:length(thesholding);
        threshold_type = thesholding(i);
        for j = 1:length(wavelets)
            wave = wavelets(j);
            for k = 1:length(p)
                rel_threhold = p(k);
                [~,~,~,threshold,SNR] = denoising_image(A_noise,wave,rel_threhold,type,threshold_type,Level,A);
                SNRMat(l,i,j,k) = SNR;
                %threshold
                ThreshMat(l,i,j,k) = threshold;
                % if (mprev > mserr)
                %     BestImage = A2;
                % end
                % mprev = mserr;
            end
            progressbar(index/n);
            index = index + 1;
        end
    end
end
save Data2.mat ThreshMat SNRMat
%% plot best image
%load Data.mat
l = 1; % 1 is redundant, 2 is standard
ibest = 1;
jbest = 1;
kbest = 1;
bestSNR = SNRMat(l,ibest,jbest,kbest)
for i = 1:size(SNRMat,2)
    for j=1:size(SNRMat,3)
        for k=1:size(SNRMat,4)
            if (bestSNR < SNRMat(l,i,j,k))
                bestSNR = SNRMat(l,i,j,k); ibest = i; jbest = j; kbest = k;
            end
        end
    end
end
BestThresh = ThreshMat(l,ibest,jbest,kbest);
[BestImage,~,~,thr,SNR] = denoising_image(A_noise,wavelets(jbest),p(kbest),TypeTransform(l),thesholding(ibest),Level,A);
% Plot the result after denosing for best params
figure
imshow(BestImage)

function [A2,mserr,mserr_rel,threshold,SNR] = denoising_image(A,wave,rel_threhold,type,threshold_type,Level,A_correct)
    % The tensor contains red, green and blue components
    Ared = A(:,:,1);
    Agreen = A(:,:,2);
    Ablue = A(:,:,3);
    
    % Recent versions of the wavelet toolbox can transform the tensor A.
    % It is generally safer to split the channels and do it manually
    %[c,l] = wavedec2(A, 4, 'bior4.4');

    [A2red,A2green,A2blue,threshold] = denoise_tensor(Ared,Agreen,Ablue,wave,rel_threhold,type,threshold_type,Level);
    %threshold
    A2 = zeros(size(A));
    A2(:,:,1) = A2red;
    A2(:,:,2) = A2green;
    A2(:,:,3) = A2blue;

    % if condition true we are denoising, thus can compute the signal to
    % noise ratio
    SNR = 10*log10(norm(A2,'fro')/(norm(A2-double(A),'fro')));
    
    A2 = uint8(A2);
    
    % Plot the result after compression
    % figure
    % imshow(A2)
    
    % Root mean square error
    %mserr = sqrt(sum(sum(sum( (A2-A_correct).^2 ))) / length(A_correct(:)));
    mserr = norm(double(A2)-double(A_correct),'fro');
    % the norm of A(:) is equal to the so-called Frobenius norm
    mserr_rel = mserr / norm(double(A_correct), "fro");
end