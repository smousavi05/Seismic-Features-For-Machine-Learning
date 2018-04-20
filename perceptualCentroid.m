function [PCmean,PCstd,PCmax,PCmin,N] = perceptualCentroid(x,fs,varargin)
% Perceptual spectral centroid
% 
%   PCMEAN = PERCEPTUALCENTROID(X,FS) calculates the spectral centroid of
%   signal X, sampled at FS, with respect to mel-frequency.
%   
%   The algorithm first calculates the spectrogram of X; each segment has
%   the maximum of length(X)/8 or 2048 samples, calculated with 50% overlap
%   and windowed with a hamming window. The spectrogram function is the
%   built-in Matlab function SPECTROGRAM. The centroid is calculated with
%   respect to mels by converting the FFT bin frequencies to mels. PCMEAN
%   is the mean of the perceptual spectral centroids calculated for each
%   segment. The perceptual spectral centroid can also be calculated using
%   ERBs, cents, or hertz (see below).
% 
%   X can be a vector, matrix, or multidimensional array;
%   PERCEPTUALCENTROID will operate along the first non-signleton
%   dimension, and return a value for each corresponding row/column/etc.
% 
%   PCmean = perceptualCentroid(x,fs,'parameter',value) allows numerous
%   parameters to be specified. These parameters are:-
%       'cref'     : {27.5} | scalar
%           Specifies the reference frequency when calulating the centroid
%           in terms of cents.
%       'dim'      : {first non-singleton dimension} | integer
%           Specify the dimension over which to calculate the perceptual
%           spectral centroid.
%       'loudness' : {'none'} | 'A' | 'B' | 'C' | 'D' | 'ISO-226'
%           Specifies whether loudness weighting is applied to the spectrum
%           prior to the centroid calculation. The default is for no
%           weighting to be applied. 'A', 'B', 'C', and 'D'  correspond to
%           frequency weighting curves defined in IEC 61672:2003; 'ISO-226'
%           applies loudness weighting derived from ISO 226:2003.
%       'nfft'     : {max([length(x)/2, 2048])} | integer
%           Specifies SPECTROGRAM's FFT size.
%       'noverlap' : {nfft/2} | integer
%           Specifies the number of samples over which SPECTROGRAM's
%           segments overlap.
%       'output'   : {'units'} | 'hz'
%           Specifies whether the output data are in Hz ('hz') or in units
%           determined by the 'scale' option ('units') (default).
%       'phon'     : {65} | scalar
%           Specifies the loudness level if using ISO 226:2003-based
%           loudness weighting.
%       'scale'    : {'mel'} | 'linear' | 'erb' | 'cents'
%           Specifies frequency scale use to calculate the centroid. 'mel'
%           uses the mel-frequency scale, 'linear' uses a linear frequency
%           scale (corresponding to the traditional spectral centroid
%           measure), 'erb' uses the equivalent-rectangular-bandwidth-rate
%           scale, and cents uses a scale based on musical cents with
%           respect to A-1 (27.5 Hz).
%       'window'   : {hamming(nfft)} | vector
%           Specifies the window applied to each SPECTROGRAM segment.
% 
%   [PCmean,PCstd,PCmax,PCmin,N] = ... returns additional information:-
%       PCmean : as described above.
%       PCstd  : the standard deviation of the perceptual spectral
%                centroids.
%       PCmax  : the maximum of the perceptual spectral centroids.
%       PCmin  : the minimum of the perceptual spectral centroids.
%       N      : the sample size for each statistic.
% 
%   For more information about the FFT calculation and its options, consult
%   the SPECTROGRAM documentation.
% 
%   Examples
% 
%     Example 1: Calculate the mel-frequency spectral centroid of random
%       noise:
% 
%       fs = 44100;
%       x = randn(fs,1);
%       PC = perceptualCentroid(x,fs);
% 
%     Example 2: Calculate the ERB-spaced spectral centroid of random
%       noise with the output in Hertz:
% 
%       PC = perceptualCentroid(x,fs,'scale','erb','output','hz');
% 
%     Example 3: Calculate the musically-spaced spectral centroid of
%       random noise using a 4096-point Hann window:
% 
%       PC = perceptualCentroid(x,fs,'scale','cents',...
%           'window',hann(4096));
% 
%   Authors: Chris Hummersone & Kirsten Hermes, 2014
%   
%   See also SPECTROGRAM.

% ==========================================================
% Last changed:     $Date: 2015-05-20 11:10:28 +0100 (Wed, 20 May 2015) $
% Last committed:   $Revision: 379 $
% Last changed by:  $Author: ch0022 $
% ==========================================================

    %% read inputs and make default assignments

    assert(~isscalar(x),'''X'' cannot be a scalar')

    dims = size(x);

    propNames = {'dim','nfft','noverlap','window','scale','output','cref','loudness','phon'}; % permitted prop names

    options = struct;
    % read parameter/value options
    if nargin>2 % if parameters are specified
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
           error('PERCEPTUALCENTROID needs propertyName/propertyValue pairs following the x and fs arguments.')
        end
        % write parameters to options struct
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
            if any(strcmpi(pair{1},propNames))
                propName = char(propNames(find(strcmpi(pair{1},propNames),1,'last')));
                options.(propName) = pair{2};
            else
                error(['Unknown option ''' pair{1} '''']);
            end
        end
    end

    % other parameters
    dim = getProperty(options,'dim',find(dims>1,1,'first'));
    assert(dim>0 && dim<=length(dims) && isscalar(dim) && round(dim)==dim,...
        '''dim'' must be greater than zero and less than or equal to the number of dimensions in X.');

    % determine window and nfft
    defWinHandle = @hamming; % default window function
    [nfft,nfft_set] = getProperty(options,'nfft',max([2048,floor(size(x,dim)/8)]));
    [window,win_set] = getProperty(options,'window',defWinHandle(nfft));

    if win_set
        % window is specified
        if isscalar(window)
            window = defWinHandle(window);
        end
        if ~nfft_set
            nfft = length(window);
        end
        if nfft~=length(window)
            error('NFFT must be equal to the window length.')
        end
    else
        % window is not specified
        window = defWinHandle(nfft);
    end

    % other parameters
    noverlap = getProperty(options,'noverlap',floor(nfft/2));
    scale = getProperty(options,'scale','mel');
    output = getProperty(options,'output','units');
    cref = getProperty(options,'cref',27.5);
    loudness = getProperty(options,'loudness','none');

    % tests
    assert(noverlap<min([nfft,length(window)]),'''noverlap'' must be less than ''nfft'' and the window length.');
    assert(isscalar(cref) && cref>0,'''cref'' must be greater than 0.');

    %% permute and rehape x to operate down columns

    order = mod(dim-1:dim+length(dims)-2,length(dims))+1;
    dims_shift = dims(order);
    x2 = rearrange(x,order,[dims_shift(1) numel(x)/dims_shift(1)]);

    %% determine functions

    fHandleDoNothing = @(x,~) (x);

    % choose frequency transformation
    switch lower(scale)
        case 'linear'
            fHandleFromF = fHandleDoNothing;
            fHandleToF = fHandleDoNothing;
        case 'mel'
            fHandleFromF = @frequency2mel;
            fHandleToF = @mel2frequency;
        case 'erb'
            fHandleFromF = @frequency2erb;
            fHandleToF = @erb2frequency;
        case 'cents'
            fHandleFromF = @frequency2cents;
            fHandleToF = @cents2frequency;
        otherwise
            error(['Requested scale ''' scale ''' not recognised. Options are ''linear'', ''mel'', ''erb'', or ''cents''.']);
    end

    % warn if cref but scale~=cents
    if ~strcmpi(scale,'cents') && ...
            any(cell2mat(strfind(lower(varargin(cellfun(@ischar,varargin))),'cref')))
        warning('Option ''cref'' only affects the output when ''scale'' is set to ''cents''.')
    end

    %% calculate centroid and stats

    % pre-allocate outputes
    PCmean = zeros(1,size(x2,2));
    PCstd = zeros(1,size(x2,2));
    PCmax = zeros(1,size(x2,2));
    PCmin = zeros(1,size(x2,2));
    N = zeros(1,size(x2,2));

    % determine output conversion
    switch lower(output)
        case 'hz'
            % use fHandleToF
        case 'units'
            fHandleToF = fHandleDoNothing;
        otherwise
            error(['Requested output ''' output ''' not recognised. Options are ''hz'', or ''units''.']);
    end

    % determine phon, if relevant
    switch lower(loudness)
        case 'iso-226'
            phon = getProperty(options,'phon',65);
        otherwise
            if isfield(options,'phon')
                warning('''phon'' option has no effect unless using ''ISO-226'' loudness weighting.')
            end
    end

    % determine loudness weighting
    switch lower(loudness)
        case 'none'
            w_l = @(x) ones(size(x));
        case 'a'
            w_l = @a_weighting;
        case 'b'
            w_l = @b_weighting;
        case 'c'
            w_l = @c_weighting;
        case 'd'
            w_l = @d_weighting;
        case 'iso-226'
            w_l = @(x) loud_weight(x,phon);
        otherwise
            error(['Requested loudness weighting ''' loudness ''' not recognised. Options are ''none'', or ''A''.']);
    end

    for c = 1:size(x2,2) % across the dim in input
        % caluclate spectrogram
        [X,f] = spectrogram(x2(:,c),window,noverlap,nfft,fs);
        % ignore frequencies greater than Nyquist
        IX = f<=fs/2;
        f = f(IX);
        % convert frequencies
        g = fHandleFromF(f,cref);
        % calculate magnitude
        mag = abs(X(IX,:));
        % calculate loudness weighting
        loud = w_l(f);
        % pre-allocate temp centroid values
        centroid = zeros(1,size(X,2));
        for d = 1:size(X,2) % iterate through spectrogram windows
            centroid(d) = fHandleToF(sum((loud.*mag(:,d)).*g)./sum(loud.*mag(:,d)),cref);
        end
        % calculate stats
        PCmean(c) = mean(centroid);
        PCstd(c) = std(centroid);
        PCmax(c) = max(centroid);
        PCmin(c) = min(centroid);
        N(c) = size(X,2);
    end

    %% inversely permute output back to input dimensions
    
    PCmean = irearrange(PCmean,order,[1 dims_shift(2:end)]);
    PCstd = irearrange(PCstd,order,[1 dims_shift(2:end)]);
    PCmax = irearrange(PCmax,order,[1 dims_shift(2:end)]);
    PCmin = irearrange(PCmin,order,[1 dims_shift(2:end)]);
    N = irearrange(N,order,[1 dims_shift(2:end)]);

end

function f = mel2frequency(m,~)
%MEL2FREQUENCY convert mel to frequency
    f = 700.*(exp(m./1127)-1);
end

function m = frequency2mel(f,~)
%FREQUENCY2MEL convert frequency to mel
    m = 1127.*log(1+(f./700));
end

function f = erb2frequency(b,~)
%ERB2FREQUENCY convert ERB to frequency
    f = (10.^(b./21.4)-1)./0.00437;
end

function b = frequency2erb(f,~)
%FREQUENCY2ERB convert frequency to ERB
    b = 21.4.*log10(1+(0.00437.*f));
end

function f = cents2frequency(c,cref)
%CENTS2FREQUENCY convert cents to frequency (ref. A-1)
    f = cref.*(2.^(c./1200)-1);
end

function c = frequency2cents(f,cref)
%FREQUENCY2CENTS convert frequency to cents (ref. A-1)
    c = 1200.*log2((f./cref)+1);
end

function [propValue,isset] = getProperty(options,propName,default)
%GETPROPERTY return property/default value
    if isfield(options,propName)
        propValue = options.(propName);
        isset = true;
    else
        propValue = default;
        isset = false;
    end
end

function w = a_weighting(f)
%A_WEIGHTING return A-weighting magnitude coefficients
    w = ((12200^2).*(f.^4))./...
        (((f.^2)+(20.6^2)).*((((f.^2)+(107.7^2)).*((f.^2)+(737.9^2))).^0.5).*((f.^2)+(12200^2)));
end

function w = b_weighting(f)
%B_WEIGHTING return B-weighting magnitude coefficients
    w = ((12200^2).*(f.^3))./...
        (((f.^2)+(20.6^2)).*sqrt((f.^2)+(158.5^2)).*((f.^2)+(12200^2)));
end

function w = c_weighting(f)
%C_WEIGHTING return C-weighting magnitude coefficients
    w = ((12200^2).*(f.^2))./...
        (((f.^2)+(20.6^2)).*((f.^2)+(12200^2)));
end

function w = d_weighting(f)
%D_WEIGHTING return D-weighting magnitude coefficients
    hf = (((1037918.48-(f.^2)).^2)+(1080768.16.*(f.^2)))./...
        (((9837328-(f.^2)).^2)+(11723776.*(f.^2)));
    w = (f./(6.8966888496476*(10^(-5)))).*sqrt(hf./(((f.^2)+79919.29).*((f.^2)+1345600)));
end

function y = rearrange(x,order,shape)
%REARRANGE reshape and permute to make target dim column
    y = permute(x,order);
    y = reshape(y,shape);
end

function y = irearrange(x,order,shape)
%IREARRANGE reshape and permute to original size
    y = reshape(x,shape);
    y = ipermute(y,order);
end
