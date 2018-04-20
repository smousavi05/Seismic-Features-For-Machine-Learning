function h = plot_sac_spectrum(fname)

% plot_sac_spectrum reads data from a sac datafile (using the function readsac)
% plots the data and the spectrum up to the Nyquist frequency
% returns a figure handle corresponding to the plot
%
% Inputs:		fname (string), name of SAC file
%
% Outputs:		h (integer) figure handle
%
%

% call readsac function with the filename
% readsac returns time and amplitude as vectors
% assign to t and a, respectively

[t,a] = readsac(fname);

% normalize seismogram amplitude and set initial time to zero

a = (a-mean(a))/max(abs(a-mean(a)));
t = t-t(1);

% plot seismogram in upper subplot

h = figure;
ax = subplot(2,1,1);
line = plot(t,a);

xlabel('Time (s)','fontname','Helvetica','fontsize',18);
ylabel('Normalized Amplitude','fontname','Helvetica','fontsize',18);
title('Normalized Seismogram','fontname','Helvetica','fontsize',18);

set(ax,'fontsize',18,'fontname','Helvetica');

set(line,'linewidth',1);

% now take spectrum

% first, set up frequency vector (goes from zero to Nyquist frequency, or half of sampling frequency)
% with length(t)/2+1 points

npoints = length(t)/2+1;

fs = 1/(t(2)-t(1));
freq = fs/2*linspace(0,1,npoints);

% take spectrum, normalizing by length of signal

spect = fft(a)/length(a);

% plot absolute value of spectrum in lower plot
% need to restrict spectrum to only npoints

ax = subplot(2,1,2);
line = loglog(freq,abs(spect(1:npoints)));

xlabel('Frequency (Hz)','fontname','Helvetica','fontsize',18);
ylabel('Spectral Amplitude','fontname','Helvetica','fontsize',18);
title('Seismogram Spectrum','fontname','Helvetica','fontsize',18);

set(ax,'fontsize',18,'fontname','Helvetica');

set(line,'linewidth',1);

end