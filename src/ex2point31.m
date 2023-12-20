% Use 1000 equispaced points on the interval [-1,1].
t = linspace(-2, 2, 1008);

% Sample a smooth function
y = abs(t) .*(2+cos(t)) .* sign(t);
% Try a non-smooth function also:
% y = abs(t) .* exp(t)
%% decontruct signal

% Compute its wavelet transform, four levels deep, using the Daubechies 2
% wavelet.
Level = 4;
wavelet = "db2";
swc = swt(y, Level, wavelet);
% Visualize the coefficients on a logarithmic scale.
figure
semilogy(abs(swc(1,:)));
%hold on
%for i = 1:Level
%    semilogy(abs(swc(i,:)));
%end
%hold off
xlabel("$i$",Interpreter="latex");
ylabel("$|c_i|$",Interpreter="latex");

%% contruct noisy signal

% apply wavelet transform on signal with added noise
rng(42)
epsilon = 5e-1;
%noise = epsilon*rand(size(y));
noise = epsilon*randn(size(y));
ynoise = y + noise;
bias = mean(noise); %mean of uniform distrubution

SNR_noise = 10*log10(norm(y,'fro')/norm(ynoise-y,'fro'));

figure
plot(t,y)
hold on
plot(t,ynoise)
hold off
xlabel("$t$",Interpreter="latex");
ylabel("function",Interpreter="latex");
legend("$f(t)$","$f(t)+\epsilon \mathcal{N}(0,1)$",Interpreter="latex")
%% denoise function
% Find small coefficents and set them to zero.
%T = max(abs(cnoise));
%delta = 10e-4*T;
p = 1e-2; type = "redundant"; thresholdType = "Soft";
wavelet = "db30";

[y2,threshold] = denoise_func(ynoise,wavelet,p,type,thresholdType,Level);
SNR = 10*log10(norm(y)/norm(y2-y));

%noiseMean = mean(abs(noise));
err = abs(y-y2);
errTotal = norm(err)
errTotalnoise = norm(noise)
%% plots

% Plot the error on a logarithmic scale. 
figure
semilogy(t, err)
hold on
%semilogy(t, noise)
yline(bias,Label="noise",Interpreter="latex")
hold off
xlabel("$t$",Interpreter="latex");
ylabel("Errors",Interpreter="latex");
legend("$|f(t_i)-\hat{f}(t_i)+mean(noise)|$",Interpreter="latex")

% plot of clean and recontructed function
figure
plot(t,y)
hold on
plot(t,y2)
hold off
xlabel("$t$",Interpreter="latex");
ylabel("function",Interpreter="latex");
legend("$f(t)$","$\hat{f}(t)-mean(noise)$",Interpreter="latex")

%% finding best threshold for certain parameters

deltaList = linspace(0,1,400);

%type = "redundant"; 
type = "other";
wavelet = "db2";
thresholdTypeList = ["Hard","Soft"];
SNRList = zeros(length(thresholdTypeList),length(deltaList));

for l = 1:length(thresholdTypeList)
    thresholdType = thresholdTypeList(l);
    for i = 1:length(deltaList)
        delta = deltaList(i);
        %[cnoise,I] = Hard_threshold(delta,cnoiseInit);
        [y2,~] = denoise_func(ynoise,wavelet,delta,type,thresholdType,Level);
        SNRList(l,i) = 10*log10(norm(y)/norm(y2-y));
    end
end

%% plots

ibest = 1;
jbest = 1;
bestSNR = SNRList(ibest,jbest);
for i = 1:size(SNRList,1)
    for j=1:size(SNRList,2)
        if (bestSNR <= SNRList(i,j))
            bestSNR = SNRList(i,j); ibest = i; jbest = j;
        end
    end
end
bestSNR

[y2, BestTr] = denoise_func(ynoise,wavelet,deltaList(jbest),type,thresholdTypeList(ibest),Level);
SNR = 10*log10(norm(y)/norm(y2-y));

% plot best recontructed function
figure
plot(t,y)
hold on
plot(t,y2)
hold off
xlabel("$t$",Interpreter="latex");
ylabel("function",Interpreter="latex");
legend("$f(t)$","$\hat{f}(t)$",Interpreter="latex")

% Plot the error on a logarithmic scale. 
err = abs(y-y2);
figure
semilogy(t, err)
%semilogy(t, noise)
hold off
xlabel("$t$",Interpreter="latex");
ylabel("Errors",Interpreter="latex");
legend("$|f(t_i)-\hat{f}(t_i)|$",Interpreter="latex")

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