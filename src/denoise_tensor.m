function [B1,B2,B3] = denoise_tensor(A1,A2,A3,wavelet,p,type,threshold_type)
    if contains(type, "redundant")
        [c1,l1,a1,b1] = swt2(A1, 4, 'bior4.4');
        [c2,l2,a2,b2] = swt2(A2, 4, 'bior4.4');
        [c3,l3,a3,b3] = swt2(A2, 4, 'bior4.4');
        T = max([max(abs(c1)), max(abs(c2)), max(abs(c3))]);
        threshold = p*T;
        if contains(threshold_type, "hard")
            [c1,I1] = Hard_threshold(threshold,c1);
            [c2,I2] = Hard_threshold(threshold,c2);
            [c3,I3] = Hard_threshold(threshold,c3);
        else
            [c1,I1] = Soft_threshold(threshold,c1);
            [c2,I2] = Soft_threshold(threshold,c2);
            [c3,I3] = Soft_threshold(threshold,c3);
        end
        B1 = iswt2(c1, l1,a1,b1, wavelet);
        B2 = iswt2(c2, l2,a2,b2, wavelet);
        B3 = iswt2(c3, l3,a3,b3, wavelet);
    else
        [c1,l1] = wavedec2(A1, 4, wavelet);
        [c2,l2] = wavedec2(A2, 4, wavelet);
        [c3,l3] = wavedec2(A3, 4, wavelet);
        T = max([max(abs(c1)), max(abs(c2)), max(abs(c3))]);
        threshold = p*T;
        if contains(threshold_type, "hard")
            [c1,I1] = Hard_threshold(threshold,c1);
            [c2,I2] = Hard_threshold(threshold,c2);
            [c3,I3] = Hard_threshold(threshold,c3);
        else
            [c1,I1] = Soft_threshold(threshold,c1);
            [c2,I2] = Soft_threshold(threshold,c2);
            [c3,I3] = Soft_threshold(threshold,c3);
        end
        
        B1 = waverec2(c1, l1, wavelet);
        B2 = waverec2(c2, l1, wavelet);
        B3 = waverec2(c3, l1, wavelet);
    end
end

function [c,I] = Hard_threshold(delta, c)
    I = find(abs(c) < delta);
    c(I) = 0;
end

function [c, I] = Soft_threshold(delta,c)
    I = find(abs(c) < delta);
    c = sign(c).*(abs(c)-delta);
    c(I) = 0;
end