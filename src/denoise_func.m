function [B,threshold] = denoise_func(A,wavelet,p,type,threshold_type,Level)
    if contains(type, "redundant",'IgnoreCase',true)
        swc = swt(A, Level, wavelet);
        T = max(abs(swc),[],"all");
        threshold = p*T;
        sz = [1,numel(swc)]; sz2 = size(swc);
        c = reshape(swc,sz);
        if contains(threshold_type, "hard",'IgnoreCase',true)
            [c,~] = Hard_threshold(threshold,c);
        else
            [c,~] = Soft_threshold(threshold,c);
        end
        swc = reshape(c,sz2);
        B = iswt(swc, wavelet);
    else
        [c,l] = wavedec(A, Level, wavelet);
        T = max(abs(c),[],"all");
        threshold = p*T;
        if contains(threshold_type, "hard",'IgnoreCase',true)
            [c,~] = Hard_threshold(threshold,c);
        else
            [c,~] = Soft_threshold(threshold,c);
        end
        B = waverec(c, l, wavelet);
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