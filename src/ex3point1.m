% The colormap is empty.
[A,cmap] = imread('bib.png');
A=im2double(A);
figure
imshow(A)

% The tensor contains red, green and blue components
Ared = A(:,:,1);
Agreen = A(:,:,2);
Ablue = A(:,:,3);
Ared_orig= Ared;


mask = zeros(size(Ared));
mask(50:60,70:80) = 1; % HARDCODED AS FUCK
mask = mask > 0;
mask_compl = ~mask;
Ared(mask) = 0;  % masked
Agreen(mask) = 0; % masked 
Ablue(mask) = 0; % masked
% Ared(mask_compl) = 0;
Ared_guess = Ared;%B0
Agreen_guess = Agreen;
Ablue_guess = Ablue;

% HARDCODED AS FUCK
Ared_extr = [Ared(49,70:80) ; Ared(61,70:80)];
Agreen_extr = [Agreen(49,70:80) ; Agreen(61,70:80)];
Ablue_extr = [Ablue(49,70:80) ; Ablue(61,70:80)];
n= size(Ared_extr,1);
Ared_interp = interp1(1:n,Ared_extr,linspace(1,n,12*(n-1)+1),'spline');
Agreen_interp = interp1(1:n,Agreen_extr,linspace(1,n,12*(n-1)+1),'spline');
Ablue_interp = interp1(1:n,Ablue_extr,linspace(1,n,12*(n-1)+1),'spline');

% HARDCODED AS FUCK
Ared_guess(50:60,70:80)=Ared_interp(2:12,:);
Agreen_guess(50:60,70:80)=Agreen_interp(2:12,:);
Ablue_guess(50:60,70:80)=Ablue_interp(2:12,:);


%% Inpainting algorithm
[Bred_prev,Bgreen_prev,Bblue_prev] = denoise_tensor(Ared_guess,Agreen_guess,Ablue_guess);
interations  = 5;
for i = 1:5
    [Bred_next,Bgreen_next,Bblue_next] = denoise_tensor(Bred_prev,Bgreen_prev,Bblue_prev);
    Bred_next(mask_compl)=0;Bgreen_next(mask_compl)=0;Bblue_next(mask_compl)=0;
    Bred_prev = Ared + Bred_next;Bgreen_prev = Agreen + Bgreen_next;Bblue_prev = Ablue + Bblue_next;
end


%% 
A2 = zeros(size(A));
A2(:,:,1) = Bred_prev;
A2(:,:,2) = Bgreen_prev;
A2(:,:,3) = Bblue_prev;
% A2 = uint8(A2);
figure
imshow(A2)