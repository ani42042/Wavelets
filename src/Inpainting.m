clear;
close all;
p=1e-2;
threshold_type="hard"
wavelets = ["db1","db4","sym4","coif4","bior4.4","rbio3.9"];
type = "redundant";
iterations = 200;
[A,cmap] = imread('bib.png');
A=imresize(A,[896 1344]);
A=im2double(A);
% figure
% imshow(A)
snr_zero =[];snr_guess =[];snr_recon = [];

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
    col_start = 649+15;col_end=659+25;
    row_start = 441+101; row_end=451+103;
    % col_start = 651;col_end=661;
    % row_start = 443; row_end=453;
    height = row_end-row_start+2;
    mask = zeros(size(Ared));
    mask(row_start:row_end,col_start:col_end) = 1;
    mask = mask > 0;
    mask_compl = ~mask;
    Ared(mask) = 0;  % masked
    Agreen(mask) = 0; % masked 
    Ablue(mask) = 0; % masked
    % Ared(mask_compl) = 0;
    Ared_guess = Ared;%B0
    Agreen_guess = Agreen;
    Ablue_guess = Ablue;
    % Ablocked = zeros(size(A));
    % Ablocked(:,:,1) = Ared;
    % Ablocked(:,:,2) = Agreen;
    % Ablocked(:,:,3) = Ablue;
    % figure
    % imshow(Ablocked);
    % [Bred_prev,Bgreen_prev,Bblue_prev] = denoise_tensor(Ared_guess,Agreen_guess,Ablue_guess,wavelet,p,type,threshold_type);
    % iterations  = 200;
    % for i = 1:iterations
    %     [Bred_next,Bgreen_next,Bblue_next] = denoise_tensor(Bred_prev,Bgreen_prev,Bblue_prev, wavelet,p,type,threshold_type);
    %     Bred_next(mask_compl)=0;Bgreen_next(mask_compl)=0;Bblue_next(mask_compl)=0;
    %     Bred_prev = Ared + Bred_next;Bgreen_prev = Agreen + Bgreen_next;Bblue_prev = Ablue + Bblue_next;
    % end
    % 'guess'
    % A2 = zeros(size(A));
    % A2(:,:,1) = Bred_prev;
    % A2(:,:,2) = Bgreen_prev;
    % A2(:,:,3) = Bblue_prev;
    % % figure
    % % imshow(A2);
    % 
    % 
    % snr_zero(w)=10*log10(norm(A,"fro")/norm(A2-A,'fro'));
    %interpolation
    Ared_extr = [Ared(row_start-1,col_start:col_end) ; Ared(row_end+1,col_start:col_end)];
    Agreen_extr = [Agreen(row_start-1,col_start:col_end) ; Agreen(row_end+1,col_start:col_end)];
    Ablue_extr = [Ablue(row_start-1,col_start:col_end) ; Ablue(row_end+1,col_start:col_end)];
    n= size(Ared_extr,1);
    Ared_interp = interp1(1:n,Ared_extr,linspace(1,n,height*(n-1)),'spline');
    Agreen_interp = interp1(1:n,Agreen_extr,linspace(1,n,height*(n-1)),'spline');
    Ablue_interp = interp1(1:n,Ablue_extr,linspace(1,n,height*(n-1)),'spline');
    
    Ared_guess(row_start:row_end,col_start:col_end)=Ared_interp(2:height,:);
    Agreen_guess(row_start:row_end,col_start:col_end)=Agreen_interp(2:height,:);
    Ablue_guess(row_start:row_end,col_start:col_end)=Ablue_interp(2:height,:);
    
    
    % Inpainting algorithm
    [Bred_prev,Bgreen_prev,Bblue_prev] = denoise_tensor(Ared_guess,Agreen_guess,Ablue_guess,wavelet,p,type,threshold_type);
    A2 = zeros(size(A));
    A2(:,:,1) = Bred_prev;
    A2(:,:,2) = Bgreen_prev;
    A2(:,:,3) = Bblue_prev;
    figure
    imshow(A2);
    snr_guess(w)=10*log10(norm(A,"fro")/norm(A2-A,'fro'));
    % iterations  = 200;
    for i = 1:iterations
    
        [Bred_next,Bgreen_next,Bblue_next] = denoise_tensor(Bred_prev,Bgreen_prev,Bblue_prev,wavelet,p,type,threshold_type);
        Bred_next(mask_compl)=0;Bgreen_next(mask_compl)=0;Bblue_next(mask_compl)=0;
        Bred_prev = Ared + Bred_next;Bgreen_prev = Agreen + Bgreen_next;Bblue_prev = Ablue + Bblue_next;
    end
    
    A2 = zeros(size(A));
    A2(:,:,1) = Bred_prev;
    A2(:,:,2) = Bgreen_prev;
    A2(:,:,3) = Bblue_prev;
    figure
    imshow(A2)
    'recon'
    snr_recon(w)=10*log10(norm(A,"fro")/norm(A2-A,'fro'));
end
% figure
% imshow(A2)
% 10*log10(norm(A,'fro')/norm(A2-A,'fro'))
%%
% Ablocked = zeros(size(A));
% Ablocked(:,:,1) = Ared;
% Ablocked(:,:,2) = Agreen;
% Ablocked(:,:,3) = Ablue;
% figure
% imshow(Ablocked);

% fid = fopen('snr_2_hard_sky.dat','wt');
% fprintf(fid,'%f ',snr_zero);
% fprintf(fid,'\n');
% fprintf(fid,'%f ',snr_guess);
% fprintf(fid,'\n');
% fprintf(fid,'%f ',snr_recon);
% fclose(fid);
