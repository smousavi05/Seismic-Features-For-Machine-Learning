clear all


% reading the waveform
[time,data,hdr] = fget_sac('master-filt');
data = data(1100:end);
time = time(1100:end);
          
% %% SFFT   
% close all
% 
%    figure(1);
%    clf;
%    set(gcf,'Position',[400,500,800,500], 'PaperPosition', [0.25 2.5 7.0 5.25]);
%    subplot 211 % Plot seismogram
%    plot(time,data);title('Bandpass filtered between 10-40 HZ')
%    max_data = max(abs(data))
%    ylim([-1.1 1.1].*max_data);
% 
%    subplot 212   
%    spectrogram(data,128,120,128,1E3,'yaxis');title('Short Time Fourier Transform');
%  % colormap winter
%    colormap(gca,1-gray);
%  % hc = colorbar('Location','SouthOutside');
% 
% %    tflog=specgram(data,256,200,32,31);
% %    imagesc(abs(tflog));


% SST yang
   figure(2);
[ss_energy coefTensor InsFreq] = SS_WP_Yang(data);
 
%% SST Dmytro Iatsenko   
% figure (3)
%[SWFT,freq,wopt] = sswft(data,100,'fmax',50,'fmin',0.01,'plot','amp++')
% title('Synchrosqueezed Windowed Fourier Transform') 


% Eugene Brevdo
   figure(3);
[Tx, fs, Wx, as, w] = synsq_cwt_fw(time, data, 64, 'morlet');
tplot(Tx, time, fs, 'Mquant')


IST=imag(Tx);
RST=real(Tx);

%% Otsu
% leveli = graythresh(IST);
% levelR = graythresh(RST);
% 
% B=fir1(5,levelR,'high');
% C=fir1(5,leveli,'high');
% 
% xfri=filter(B,1,RST)+j*filter(C,1,IST); 
% X = synsq_cwt_iw(xfri, 64, 'morlet');
% plot(time,X);

% function STFTcoef = STFT(f, time_win, factor_redund, f_sampling)
% %
% % 1D Windowed Fourier Transform. 
% %
% % Input:
% % - f: Input 1D signal.
% % - time_win: window size in time (in millisecond).
% % - factor_redund: logarithmic redundancy factor. The actual redundancy
% %   factor is 2^factor_redund. When factor_redund=1, it is the minimum
% %   twice redundancy. 
% % - f_sampling: the signal sampling frequency in Hz.
% %
% % Output:
% % - STFTcoef: Spectrogram. Column: frequency axis from -pi to pi. Row: time
% %   axis. 
% %
% % Remarks:
% % 1. The last few samples at the end of the signals that do not compose a complete
% %    window are ignored in the transform in this Version 1. 
% % 2. Note that the reconstruction will not be exact at the beginning and
% %    the end of, each of half window size. However, the reconstructed
% %    signal will be of the same length as the original signal. 
% %
% % See also:
% % inverseSTFT
% %
% % Guoshen Yu
% % Version 1, Sept 15, 2006
% 
% 
% % Check that f is 1D
% if length(size(f)) ~= 2 | (size(f,1)~=1 && size(f,2)~=1)
%     error('The input signal must 1D.');
% end
% 
% if size(f,2) == 1
%     f = f';
% end
% 
% % Window size
% size_win = round(time_win/1000 * f_sampling);
% 
% % Odd size for MakeHanning
% if mod(size_win, 2) == 0
%     size_win = size_win + 1;
% end
% halfsize_win =  (size_win - 1) / 2;
% 
% w_hanning = MakeHanning(size_win); 
% 
% Nb_win = floor(length(f) / size_win * 2);
% 
% % STFTcoef = zeros(2^(factor_redund-1), size_win, Nb_win-1);
% STFTcoef = zeros(size_win, (2^(factor_redund-1) * Nb_win-1));
% 
% shift_k = round(halfsize_win / 2^(factor_redund-1));
% % Loop over 
% for k = 1 : 2^(factor_redund-1)    
%     % Loop over windows
%     for j = 1 : Nb_win - 2 % Ingore the last few coefficients that do not make a window
%         f_win = f(shift_k*(k-1)+(j-1)*halfsize_win+1 : shift_k*(k-1)+(j-1)*halfsize_win+size_win);
%         STFTcoef(:, (k-1)+2^(factor_redund-1)*j) = fft(f_win .* w_hanning');
%     end
% end
% 
% 
% function f_rec = inverseSTFT(STFTcoef, time_win, factor_redund, f_sampling, length_f)
% %
% % Inverse windowed Fourier transform. 
% %
% % Input:
% % - STFTcoef: Spectrogram. Column: frequency axis from -pi to pi. Row: time
% %   axis. (Output of STFT). 
% % - time_win: window size in time (in millisecond).
% % - factor_redund: logarithmic redundancy factor. The actual redundancy
% %   factor is 2^factor_redund. When factor_redund=1, it is the minimum
% %   twice redundancy. 
% % - f_sampling: the signal sampling frequency in Hz.
% % - length_f: length of the signal. 
% %
% % Output:
% % - f_rec: reconstructed signal. 
% %
% % Remarks:
% % 1. The last few samples at the end of the signals that do not compose a complete
% %    window are ignored in the forward transform STFT of Version 1. 
% % 2. Note that the reconstruction will not be exact at the beginning and
% %    the end of, each of half window size. 
% %
% % See also:
% % STFT
% %
% % Guoshen Yu
% % Version 1, Sept 15, 2006
% 
% % Window size
% size_win = round(time_win/1000 * f_sampling);
% 
% % Odd size for MakeHanning
% if mod(size_win, 2) == 0
%     size_win = size_win + 1;
% end
% halfsize_win =  (size_win - 1) / 2;
% 
% Nb_win = floor(length_f / size_win * 2);
% 
% % Reconstruction
% f_rec = zeros(1, length_f);
% 
% shift_k = round(halfsize_win / 2^(factor_redund-1));
% 
% % Loop over windows 
% for k = 1 : 2^(factor_redund-1)
%     for j = 1 : Nb_win - 1
%         f_win_rec = ifft(STFTcoef(:, (k-1)+2^(factor_redund-1)*j));
%         f_rec(shift_k*(k-1)+(j-1)*halfsize_win+1 : shift_k*(k-1)+(j-1)*halfsize_win+size_win) =  f_rec(shift_k*(k-1)+(j-1)*halfsize_win+1 : shift_k*(k-1)+(j-1)*halfsize_win+size_win) +  (f_win_rec');
%     end
% end
% 
% f_rec = f_rec / 2^(factor_redund-1);

%    subplot 313 % FFT plot
%   NFFT = 2^nextpow2(length(data));
%   Y=fft(sig,NFFT)/length(data);
%   FS=1/hdr.times.delta;
%   f=FS/2*linspace(0,1,NFFT/2+1)
% 
%   plot(f,2*abs(data(1:NFFT/2+1)));
%   title('Amplitude Spectrum')
%   xlabel('Frequency (Hz)')
%   ylabel('abs signal')