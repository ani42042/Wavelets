% Read the image. In this case, this results in a
% M x N x 3 tensor, with RGB components.
% The colormap is empty.
[A,cmap] = imread('bib.png');
figure
imshow(A)

%% compress image
p = 1e-3; thresholdType = 'Hard'; waveletType = 'bior4.4';
A_noise = A;
[A2,mserr,mserr_rel,compression_ratio,threshold,SNR] = denoising_image(A_noise,p,thresholdType,waveletType,A);

% Plot the result after denosing/compression
figure
imshow(A2)

%% Add noise to image

% add noise to image
rng(42)

%A_noise = imnoise(A,'salt & pepper');
%A_noise = imnoise(A,'speckle');
A_noise = imnoise(A,"gaussian");

msNoise = sqrt(sum(sum(sum( (A_noise-A).^2 ))) / length(A(:)));

% Plot the noisy image
figure
imshow(A_noise)

%% Denoise image
p = 5e-3; thresholdType = 'Soft'; waveletType = 'db1';%'bior4.4'; 
[A2,mserr,mserr_rel,compression_ratio,threshold,SNR] = denoising_image(A_noise,p,thresholdType,waveletType,A);

% Plot the result after denosing/compression
figure
imshow(A2)

function [c,I] = Hard_threshold(delta, c)
    I = find(abs(c) < delta);
    c(I) = 0;
end

function [c, I] = Soft_threshold(delta,c)
    I = find(abs(c) < delta);
    c = sign(c).*(abs(c)-delta);
    c(I) = 0;
end

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
    %I1 = find(abs(c1) < threshold);
    %c1(I1) = 0;
    %I2 = find(abs(c2) < threshold);
    %c2(I2) = 0;
    %I3 = find(abs(c3) < threshold);
    %c3(I3) = 0;
    if (strcmp('Hard',threshtype))
        [c1,I1] = Hard_threshold(threshold, c1);
        [c2,I2] = Hard_threshold(threshold, c2);
        [c3,I3] = Hard_threshold(threshold, c3);
    else
        [c1,I1] = Soft_threshold(threshold, c1);
        [c2,I2] = Soft_threshold(threshold, c2);
        [c3,I3] = Soft_threshold(threshold, c3);
    end
    A2red = waverec2(c1, l1, waveletType);
    A2green = waverec2(c2, l1, waveletType);
    A2blue = waverec2(c3, l1, waveletType);
    A2 = zeros(size(A));
    A2(:,:,1) = A2red;
    A2(:,:,2) = A2green;
    A2(:,:,3) = A2blue;

    % if condition true we are denoising, thus can compute the signal to
    % noise ratio
    if nargin >= 5
        SNR = 10*log10(norm(A2,'fro')/(norm(A2-double(A),'fro')));
    end
    
    A2 = uint8(A2);
    
    % Plot the result after compression
    % figure
    % imshow(A2)
    
    % Compression ratio
    N = length(c1)+length(c2)+length(c3);
    M = length(I1)+length(I2)+length(I3);
    compression_ratio = N / (N-M);
    
    % Root mean square error
    mserr = sqrt(sum(sum(sum( (A2-A_correct).^2 ))) / length(A_correct(:)));
    % the norm of A(:) is equal to the so-called Frobenius norm
    mserr_rel = mserr / norm(double(A_correct), "fro");

end