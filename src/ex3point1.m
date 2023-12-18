clear;
close all;
p=1e-2;
wavelets = ["db1","db4","sym4", "bior4.4","rbio3.9"];
% wavelet = wavelets(4)
% The colormap is empty.
[A,cmap] = imread('bib.png');
A=im2double(A);
% figure
% imshow(A)
snr_guess =[];snr_recon = [];
% The tensor contains red, green and blue components
%%
for w = 1:length(wavelets) 
    wavelet = wavelets(w)
    Ared = A(:,:,1);
    Agreen = A(:,:,2);
    Ablue = A(:,:,3);
    Ared_orig= Ared;
    
    % Preprocessing 
    % col_start = 1300;col_end=1310;
    % row_start = 50; row_end=60;
    col_start = 649+16;col_end=659+24;
    row_start = 441+100; row_end=451+102;
    % col_start = 550;col_end=659+150;
    % row_start = 441+100-6; row_end=451+102+5;
    height = row_end-row_start+2;
    mask = zeros(size(Ared));
    mask(row_start:row_end,col_start:col_end) = 1; % SOFTCODED AS FUCK
    mask = mask > 0;
    mask_compl = ~mask;
    Ared(mask) = 0;  % masked
    Agreen(mask) = 0; % masked 
    Ablue(mask) = 0; % masked
    % Ared(mask_compl) = 0;
    Ared_guess = Ared;%B0
    Agreen_guess = Agreen;
    Ablue_guess = Ablue;
    
    [Bred_prev,Bgreen_prev,Bblue_prev] = denoise_tensor(Ared_guess,Agreen_guess,Ablue_guess,wavelet,p);
    iterations  = 50;
    for i = 1:iterations
        [Bred_next,Bgreen_next,Bblue_next] = denoise_tensor(Bred_prev,Bgreen_prev,Bblue_prev, wavelet,p);
        Bred_next(mask_compl)=0;Bgreen_next(mask_compl)=0;Bblue_next(mask_compl)=0;
        Bred_prev = Ared + Bred_next;Bgreen_prev = Agreen + Bgreen_next;Bblue_prev = Ablue + Bblue_next;
    end
    'guess'
    A2 = zeros(size(A));
    A2(:,:,1) = Bred_prev;
    A2(:,:,2) = Bgreen_prev;
    A2(:,:,3) = Bblue_prev;
    snr_guess(w)=10*log10(norm(A,"fro")/norm(A2-A,'fro'));
    % SOFTCODED AS FUCK
    Ared_extr = [Ared(row_start-1,col_start:col_end) ; Ared(row_end+1,col_start:col_end)];
    Agreen_extr = [Agreen(row_start-1,col_start:col_end) ; Agreen(row_end+1,col_start:col_end)];
    Ablue_extr = [Ablue(row_start-1,col_start:col_end) ; Ablue(row_end+1,col_start:col_end)];
    n= size(Ared_extr,1);
    Ared_interp = interp1(1:n,Ared_extr,linspace(1,n,height*(n-1)),'spline');
    Agreen_interp = interp1(1:n,Agreen_extr,linspace(1,n,height*(n-1)),'spline');
    Ablue_interp = interp1(1:n,Ablue_extr,linspace(1,n,height*(n-1)),'spline');
    
    % SOFTCODED AS FUCK
    Ared_guess(row_start:row_end,col_start:col_end)=Ared_interp(2:height,:);
    Agreen_guess(row_start:row_end,col_start:col_end)=Agreen_interp(2:height,:);
    Ablue_guess(row_start:row_end,col_start:col_end)=Ablue_interp(2:height,:);
    
    
    % Inpainting algorithm
    [Bred_prev,Bgreen_prev,Bblue_prev] = denoise_tensor(Ared_guess,Agreen_guess,Ablue_guess,wavelet,p);
    iterations  = 50;
    for i = 1:iterations
    
        [Bred_next,Bgreen_next,Bblue_next] = denoise_tensor(Bred_prev,Bgreen_prev,Bblue_prev,wavelet,p);
        Bred_next(mask_compl)=0;Bgreen_next(mask_compl)=0;Bblue_next(mask_compl)=0;
        Bred_prev = Ared + Bred_next;Bgreen_prev = Agreen + Bgreen_next;Bblue_prev = Ablue + Bblue_next;
    end
    
    A2 = zeros(size(A));
    A2(:,:,1) = Bred_prev;
    A2(:,:,2) = Bgreen_prev;
    A2(:,:,3) = Bblue_prev;
    % A2 = uint8(A2);
    'recon'
    snr_recon(w)=10*log10(norm(A,"fro")/norm(A2-A,'fro'));
end
% figure
% imshow(A2)
% 10*log10(norm(A,'fro')/norm(A2-A,'fro'))
%%
Ablocked = zeros(size(A));
Ablocked(:,:,1) = Ared;
Ablocked(:,:,2) = Agreen;
Ablocked(:,:,3) = Ablue;
figure
imshow(Ablocked);
