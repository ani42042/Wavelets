function [B1,B2,B3,threshold] = denoise_tensor(A1,A2,A3,wavelet,p,type,threshold_type,Level)
    if nargin < 8
        Level = 4;
    end

    if contains(type, "redundant",'IgnoreCase',true)
        swc1 = swt2(A1, Level, wavelet);
        swc2 = swt2(A2, Level, wavelet);
        swc3 = swt2(A3, Level, wavelet);
        %c1 = c1T(:,:,3); c2 = c2T(:,:,3); c3 = c3T(:,:,3);
        T = max([max(abs(swc1),[],"all"), max(abs(swc2),[],'all'), max(abs(swc3),[],'all')]);
        threshold = p*T;
        sz = [1,numel(swc1)]; sz2 = size(swc1);
        c1 = reshape(swc1,sz); c2 = reshape(swc2,sz); c3 = reshape(swc3,sz);
        if contains(threshold_type, "hard",'IgnoreCase',true)
            [c1,~] = Hard_threshold(threshold,c1);
            [c2,~] = Hard_threshold(threshold,c2);
            [c3,~] = Hard_threshold(threshold,c3);
        else
            [c1,~] = Soft_threshold(threshold,c1);
            [c2,~] = Soft_threshold(threshold,c2);
            [c3,~] = Soft_threshold(threshold,c3);
        end
        swc1 = reshape(c1,sz2); swc2 = reshape(c2,sz2); swc3 = reshape(c3,sz2);
        %c1T(:,:,1) = c1; c2T(:,:,1) = c2; c3T(:,:,1) = c3;
        zerodetails = zeros(size(c1));
        B1 = iswt2(swc1, wavelet);
        B2 = iswt2(swc2, wavelet);
        B3 = iswt2(swc3, wavelet);
    else
        [c1,l1] = wavedec2(A1, Level, wavelet);
        [c2,l2] = wavedec2(A2, Level, wavelet);
        [c3,l3] = wavedec2(A3, Level, wavelet);
        T = max([max(abs(c1),[],"all"), max(abs(c2),[],'all'), max(abs(c3),[],'all')]);
        threshold = p*T;
        if contains(threshold_type, "hard",'IgnoreCase',true)
            [c1,~] = Hard_threshold(threshold,c1);
            [c2,~] = Hard_threshold(threshold,c2);
            [c3,~] = Hard_threshold(threshold,c3);
        else
            [c1,~] = Soft_threshold(threshold,c1);
            [c2,~] = Soft_threshold(threshold,c2);
            [c3,~] = Soft_threshold(threshold,c3);
        end
        
        B1 = waverec2(c1, l1, wavelet);
        B2 = waverec2(c2, l2, wavelet);
        B3 = waverec2(c3, l3, wavelet);
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