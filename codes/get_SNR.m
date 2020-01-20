function SNR = get_SNR(fo, f)

if (sum(fo(:).^2) == 0) | (sum((fo(:)-f(:)).^2) == 0)
    SNR = 0;
else
    SNR = 10*log10(sum(fo(:).^2) / (sum((fo(:)-f(:)).^2)));
end