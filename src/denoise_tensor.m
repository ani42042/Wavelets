function [B1,B2,B3] = denoise_tensor(A1,A2,A3,wavelet,p)
    [c1,l1] = wavedec2(A1, 4, wavelet);
    [c2,l2] = wavedec2(A2, 4, wavelet);
    [c3,l3] = wavedec2(A3, 4, wavelet);
    
    % Compress by putting small coefficients to zero.
    % Experiment with the threshold p.
    T = max([max(abs(c1)), max(abs(c2)), max(abs(c3))]);
    threshold = p*T;
    I1 = find(abs(c1) < threshold);
    c1(I1) = 0;
    I2 = find(abs(c2) < threshold);
    c2(I2) = 0;
    I3 = find(abs(c3) < threshold);
    c3(I3) = 0;
    
    B1 = waverec2(c1, l1, wavelet);
    B2 = waverec2(c2, l1, wavelet);
    B3 = waverec2(c3, l1, wavelet);
end