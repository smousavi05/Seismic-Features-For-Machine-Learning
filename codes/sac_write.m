function sac_write( filename, waveform, varargin )
%SAC_WRITE( filename, waveform, varargin )
%
%  This function writes waveform data and headers to a SAC file.
%
%  Example:
%      SAC_WRITE('file.sac', waveform, 'stla', 35.1232, 'stlo', -89.9320, 'idep', 'iacc', 'cmpaz', 'Z', cmpinc, 'Z');
%      SAC_WRITE('file.sac', waveform, 'idep', 'ivel', 'cmpaz', 'E', cmpinc, 'E');
%      SAC_WRITE('file.sac', waveform, 'idep', 6, 'cmpaz', 90, cmpinc, 90);
%
%   NOTE: the following variables -must- be specified to create a SAC file
%     DELTA    the sampling rate of the data in seconds
%     IDEP     the type of data (e.g., displacement, velocity, acceleration..
%                see IRIS website for format)
%
%     Acceptable values for IDEP:
%       - IUNKN  (Unknown)
%       - IDISP  (Displacement in nm)
%       - IVEL   (Velocity in nm/sec)
%       - IVOLTS (Velocity in volts)
%       - IACC   (Acceleration in nm/sec/sec)
%
%   Additional information:
%     - Many fields may be specified by either their string or their integer
%       value. For more information, see LINK 1, below.
%
%     - CMPAZ and CMPINC may be specified by either single precision
%       floating point degrees or by 'Z', 'N', or 'E'. SAC_WRITE will
%       automatically convert 'Z', 'N', and 'E' into the appropriate degrees.
%
%
%   More information on the SAC file format can be found at:
%     1. http://www.iris.edu/manuals/sac/SAC_Manuals/FileFormatPt2.html
%     2. http://www.iris.edu/software/sac/manual/file_format.html

% Version: 1.31
% Date: 2012/09/26 11:43:22
% Author(s): Brian Young

% CHANGELOG
%   1.31
%     - string headers are now properly padded with spaces

% ----------------------------------------------------------------------- %

%{
% make sure the output file doesn't already exist
if exist(filename,'file')
    reply = input(sprintf('File %s already exists. Overwrite? [y/n]: ',filename), 's');
    if (lower(reply) ~= 'y')
        disp('sac_write aborted.');
        return;
    end
end
%}

% make sure that our waveform data is an Nx1 or 1xN matrix
assert( isvector(waveform), 'ERROR: ''waveform'' should be an Nx1 or 1xN matrix.');

% make sure there are an even number of input vars
assert( ~mod( numel(varargin), 2 ) , ...
    'Error in the number of input variables. Did you forget to set a header?');

% these are the required variables
req_vars = {'DELTA', 'IDEP'};

% check to see if the required variables are in the input arguments
check_vars = ~ismember( lower(req_vars), lower(varargin(1:2:end)) );

% make sure the required variables are in the input arguments, or give error
assert( ~any(check_vars), ...
    [ sprintf('Some required variables have not been specified:') ...
    sprintf(' %s',req_vars{check_vars}) ...
    sprintf('\nSee ''doc sac_write'' for more information.')] );


% --------------------------------------------------------------- %
%        INFORMATION about the SAC headers and data format        %
%     http://www.iris.edu/software/sac/manual/file_format.html    %
%  http://www.iris.edu/manuals/sac/SAC_Manuals/FileFormatPt2.html %
% --------------------------------------------------------------- %

% pre-allocate header
h = struct( ...
    'float32', single( repmat(-12345.0, 14, 5) ), ...  % single precision
    'int32', int32( repmat(-12345, 8, 5) ), ...  % integer
    'char', char( repmat('-12345..', 8, 3) ) );  % string
h.int32(8,:) = false(1,5);  % logical (i.e., true = 1, false = 0)


% set defaults
% i've comment out the header values that by default are blank (i.e., -12345)
% h.float32( 1,1) = -12345.0;  % delta      time increment
h.float32( 1,2) = min(waveform);   % depmin     minimum amplitude
h.float32( 1,3) = max(waveform);   % depmax     maximum amplitude
% h.float32( 1,4) = -12345.0;  % scale      multiplying scale factor [not currently used]
% h.float32( 1,5) = -12345.0;  % odelta     observed time increment
h.float32( 2,1) = 0;           % b          begin time
% h.float32( 2,2) = -12345.0;  % e          end time
% h.float32( 2,3) = -12345.0;  % o          event origin marker
% h.float32( 2,4) = -12345.0;  % a          first arrival (P) marker
% h.float32( 2,5) = -12345.0;  % internal1  SAC internal variable
% h.float32( 3,1) = -12345.0;  % t0         time pick 0 (S) marker
% h.float32( 3,2) = -12345.0;  % t1         user-defined time pick 1
% h.float32( 3,3) = -12345.0;  % t2         user-defined time pick 2
% h.float32( 3,4) = -12345.0;  % t3         user-defined time pick 3
% h.float32( 3,5) = -12345.0;  % t4         user-defined time pick 4
% h.float32( 4,1) = -12345.0;  % t5         user-defined time pick 5
% h.float32( 4,2) = -12345.0;  % t6         user-defined time pick 6
% h.float32( 4,3) = -12345.0;  % t7         user-defined time pick 7
% h.float32( 4,4) = -12345.0;  % t8         user-defined time pick 8
% h.float32( 4,5) = -12345.0;  % t9         user-defined time pick 9
% h.float32( 5,1) = -12345.0;  % f          end of event time
% h.float32( 5,2) = -12345.0;  % resp0      intrument response parameter 0 [not currently used]
% h.float32( 5,3) = -12345.0;  % resp1      intrument response parameter 1 [not currently used]
% h.float32( 5,4) = -12345.0;  % resp2      intrument response parameter 2 [not currently used]
% h.float32( 5,5) = -12345.0;  % resp3      intrument response parameter 3 [not currently used]
% h.float32( 6,1) = -12345.0;  % resp4      intrument response parameter 4 [not currently used]
% h.float32( 6,2) = -12345.0;  % resp5      intrument response parameter 5 [not currently used]
% h.float32( 6,3) = -12345.0;  % resp6      intrument response parameter 6 [not currently used]
% h.float32( 6,4) = -12345.0;  % resp7      intrument response parameter 7 [not currently used]
% h.float32( 6,5) = -12345.0;  % resp8      intrument response parameter 8 [not currently used]
% h.float32( 7,1) = -12345.0;  % resp9      intrument response parameter 9 [not currently used]
% h.float32( 7,2) = -12345.0;  % stla       station latitude
% h.float32( 7,3) = -12345.0;  % stlo       station longitude
% h.float32( 7,4) = -12345.0;  % stel       station elevation
% h.float32( 7,5) = -12345.0;  % stdp       station depth [not currently used]
% h.float32( 8,1) = -12345.0;  % evla       event latitude
% h.float32( 8,2) = -12345.0;  % evlo       event longitude
% h.float32( 8,3) = -12345.0;  % evel       event elevation [not currently used]
% h.float32( 8,4) = -12345.0;  % evdp       event depth
% h.float32( 8,5) = -12345.0;  % mag        event magnitude
% h.float32( 9,1) = -12345.0;  % user0      user-defined variable 0
% h.float32( 9,2) = -12345.0;  % user1      user-defined variable 1
% h.float32( 9,3) = -12345.0;  % user2      user-defined variable 2
% h.float32( 9,4) = -12345.0;  % user3      user-defined variable 3
% h.float32( 9,5) = -12345.0;  % user4      user-defined variable 4
% h.float32(10,1) = -12345.0;  % user5      user-defined variable 5
% h.float32(10,2) = -12345.0;  % user6      user-defined variable 6
% h.float32(10,3) = -12345.0;  % user7      user-defined variable 7
% h.float32(10,4) = -12345.0;  % user8      user-defined variable 8
% h.float32(10,5) = -12345.0;  % user9      user-defined variable 9
% h.float32(11,1) = -12345.0;  % dist       source receiver distance (km)
% h.float32(11,2) = -12345.0;  % az         azimuth
% h.float32(11,3) = -12345.0;  % baz        back azimuth
% h.float32(11,4) = -12345.0;  % gcarc      great circle distance (deg)
% h.float32(11,5) = -12345.0;  % internal2  SAC internal variable
% h.float32(12,1) = -12345.0;  % internal3  SAC internal variable
h.float32(12,2) = mean(waveform);  % depmen     mean amplitude
% h.float32(12,3) = -12345.0;  % cmpaz      component azimuth
% h.float32(12,4) = -12345.0;  % cmpinc     component incident angle
% h.float32(12,5) = -12345.0;  % xminimum   minimum value of x (spectral files only)
% h.float32(13,1) = -12345.0;  % xmaximum   maximum value of x (spectral files only)
% h.float32(13,2) = -12345.0;  % yminimum   minimum value of y (spectral files only)
% h.float32(13,3) = -12345.0;  % ymaximum   maximum value of y (spectral files only)
% h.float32(13,4) = -12345.0;  % unused1    SAC unused header space
% h.float32(13,5) = -12345.0;  % unused2    SAC unused header space
% h.float32(14,1) = -12345.0;  % unused3    SAC unused header space
% h.float32(14,2) = -12345.0;  % unused4    SAC unused header space
% h.float32(14,3) = -12345.0;  % unused5    SAC unused header space
% h.float32(14,4) = -12345.0;  % unused6    SAC unused header space
% h.float32(14,5) = -12345.0;  % unused7    SAC unused header space

% h.int32(1,1) = -12345;     % nzyear     GMT year
% h.int32(1,2) = -12345;     % nzjday     GMT julian day
% h.int32(1,3) = -12345;     % nzhour     GMT hour
% h.int32(1,4) = -12345;     % nzmin      GMT minute
% h.int32(1,5) = -12345;     % nzsec      GMT second
% h.int32(2,1) = -12345;     % nzmsec     GMT millisecond
h.int32(2,2) = 6;            % nvhdr      header version number
% h.int32(2,3) = -12345;     % norid      origin id
% h.int32(2,4) = -12345;     % nevid      event id
h.int32(2,5) = numel(waveform);  % npts       number of data points
% h.int32(3,1) = -12345;     % internal4  SAC internal variable
% h.int32(3,2) = -12345;     % nwfid      waveform id
% h.int32(3,3) = -12345;     % nxsize     spectral length (spectral files only)
% h.int32(3,4) = -12345;     % nysize     spectral width (spectral files only)
% h.int32(3,5) = -12345;     % unused8    SAC unused header space
h.int32(4,1) = sac_ascii_value('itime');            % iftype     type of file
% h.int32(4,2) = -12345;     % idep       type of independent variable
h.int32(4,3) = sac_ascii_value('ib');            % iztype     reference time equivalence
% h.int32(4,4) = -12345;     % unused9    SAC unused header space
% h.int32(4,5) = -12345;     % iinst      type of recording instrument [not currently used]
% h.int32(5,1) = -12345;     % istreg     station geographic region [not currently used]
% h.int32(5,2) = -12345;     % ievreg     event geographic region [not currently used]
% h.int32(5,3) = -12345;     % ievtyp     type of event
% h.int32(5,4) = -12345;     % iqual      quality of data [not currently used]
% h.int32(5,5) = -12345;     % isynth     synthetic data flag [not currently used]
% h.int32(6,1) = -12345;     % imagtyp    magnitude type
% h.int32(6,2) = -12345;     % imagsrc    source of magnitude information
% h.int32(6,3) = -12345;     % unused10   SAC unused header space
% h.int32(6,4) = -12345;     % unused11   SAC unused header space
% h.int32(6,5) = -12345;     % unused12   SAC unused header space
% h.int32(7,1) = -12345;     % unused13   SAC unused header space
% h.int32(7,2) = -12345;     % unused14   SAC unused header space
% h.int32(7,3) = -12345;     % unused15   SAC unused header space
% h.int32(7,4) = -12345;     % unused16   SAC unused header space
% h.int32(7,5) = -12345;     % unused17   SAC unused header space
h.int32(8,1) = 1;            % leven      TRUE if waveform data is evenly spaced
h.int32(8,2) = 1;            % lpspol     TRUE if station components have + polarity
h.int32(8,3) = 1;            % lovrok     TRUE if it's OK to overwrite this file
h.int32(8,4) = 1;            % lcalda     TRUE if DIST, AZ, BAZ, and GCARC are to be auto-calculated from station and event coordinates.
% h.int32(8,5) = -12345;     % unused18   SAC unused header space

% h.char(1, 1: 8) = '-12345..';          % kstnm    station name
% h.char(1, 9:24) = '-12345..-12345..';  % kevnm    event name
% h.char(2, 1: 8) = '-12345..';          % khole    hole identification if nuclear event
h.char(2, 9:16) = 'origin  ';            % ko       event origin time string
% h.char(2,17:24) = '-12345..';          % ka       first arrival time string
% h.char(3, 1: 8) = '-12345..';          % kt0      user-defined pick string 0
% h.char(3, 9:16) = '-12345..';          % kt1      user-defined pick string 1
% h.char(3,17:24) = '-12345..';          % kt2      user-defined pick string 2
% h.char(4, 1: 8) = '-12345..';          % kt3      user-defined pick string 3
% h.char(4, 9:16) = '-12345..';          % kt4      user-defined pick string 4
% h.char(4,17:24) = '-12345..';          % kt5      user-defined pick string 5
% h.char(5, 1: 8) = '-12345..';          % kt6      user-defined pick string 6
% h.char(5, 9:16) = '-12345..';          % kt7      user-defined pick string 7
% h.char(5,17:24) = '-12345..';          % kt8      user-defined pick string 8
% h.char(6, 1: 8) = '-12345..';          % kt9      user-defined pick string 9
% h.char(6, 9:16) = '-12345..';          % kf       end of event time string
% h.char(6,17:24) = '-12345..';          % kuser0   user-defined variable string 0
% h.char(7, 1: 8) = '-12345..';          % kuser1   user-defined variable string 1
% h.char(7, 9:16) = '-12345..';          % kuser2   user-defined variable string 2
% h.char(7,17:24) = '-12345..';          % kcmpnm   component name
% h.char(8, 1: 8) = '-12345..';          % knetwk   name of seismic network
% h.char(8, 9:16) = '-12345..';          % datrd    date data was read onto computer
% h.char(8,17:24) = '-12345..';          % kinst    generic name of instrument

% get all the input arguments
for n = 1 : 2 : numel(varargin)
    
    tmp.hdr = lower( varargin{n} );    % header name (e.g., 'stla')
    tmp.val = varargin{n+1};  % header value (e.g., 35.12)
    
    switch tmp.hdr
        % sometimes we might just want to pass an entire header block as a
        % variable instead of each header individually
        % e.g., SAC_WRITE('file.sac', waveform, 'float32', float_hdr, ...)
        case 'float32'    % take an entire header block has an input (e.g., 'float32' from SAC_READ)
            assertions('float32',tmp.hdr,tmp.val);
            h.float32 = tmp.val;
            
        case 'int32'    % take an entire header block has an input (e.g., 'int32' from SAC_READ)
            assertions('int32',tmp.hdr,tmp.val);
            h.int32 = tmp.val;
            
        case 'char'    % take an entire header block has an input (e.g., 'char' from SAC_READ)
            assertions('char',tmp.hdr,tmp.val);
            h.char = tmp.val;
            
        case 'delta'    % time increment
            assertions('number',tmp.hdr,tmp.val);
            h.float32(1,1) = tmp.val;
            
        case 'depmin'    % minimum amplitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(1,2) = tmp.val;
            
        case 'depmax'    % maximum amplitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(1,3) = tmp.val;
            
            %{
        case 'scale'    % multiplying scale factor [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(1,4) = tmp.val;
            %}
            
        case 'odelta'    % observed time increment
            assertions('number',tmp.hdr,tmp.val);
            h.float32(1,5) = tmp.val;
            
        case 'b'    % begin time
            assertions('number',tmp.hdr,tmp.val);
            h.float32(2,1) = tmp.val;
            
        case 'e'    % end time
            assertions('number',tmp.hdr,tmp.val);
            h.float32(2,2) = tmp.val;
            
        case 'o'    % event origin marker
            assertions('number',tmp.hdr,tmp.val);
            h.float32(2,3) = tmp.val;
            
        case 'a'    % first arrival (P) marker
            assertions('number',tmp.hdr,tmp.val);
            h.float32(2,4) = tmp.val;
            
            %{
        case 'internal1'    % SAC internal variable
            assertions('number',tmp.hdr,tmp.val);
            h.float32(2,5) = tmp.val;
            %}
            
        case 't0'    % time pick 0 (S) marker
            assertions('number',tmp.hdr,tmp.val);
            h.float32(3,1) = tmp.val;
            
        case 't1'    % user-defined time pick 1
            assertions('number',tmp.hdr,tmp.val);
            h.float32(3,2) = tmp.val;
            
        case 't2'    % user-defined time pick 2
            assertions('number',tmp.hdr,tmp.val);
            h.float32(3,3) = tmp.val;
            
        case 't3'    % user-defined time pick 3
            assertions('number',tmp.hdr,tmp.val);
            h.float32(3,4) = tmp.val;
            
        case 't4'    % user-defined time pick 4
            assertions('number',tmp.hdr,tmp.val);
            h.float32(3,5) = tmp.val;
            
        case 't5'    % user-defined time pick 5
            assertions('number',tmp.hdr,tmp.val);
            h.float32(4,1) = tmp.val;
            
        case 't6'    % user-defined time pick 6
            assertions('number',tmp.hdr,tmp.val);
            h.float32(4,2) = tmp.val;
            
        case 't7'    % user-defined time pick 7
            assertions('number',tmp.hdr,tmp.val);
            h.float32(4,3) = tmp.val;
            
        case 't8'    % user-defined time pick 8
            assertions('number',tmp.hdr,tmp.val);
            h.float32(4,4) = tmp.val;
            
        case 't9'    % user-defined time pick 9
            assertions('number',tmp.hdr,tmp.val);
            h.float32(4,5) = tmp.val;
            
        case 'f'    % end of event time
            assertions('number',tmp.hdr,tmp.val);
            h.float32(5,1) = tmp.val;
            
            %{
        case 'resp0'    % intrument response parameter 0 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(5,2) = tmp.val;
            
        case 'resp1'    % intrument response parameter 1 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(5,3) = tmp.val;
            
        case 'resp2'    % intrument response parameter 2 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(5,4) = tmp.val;
            
        case 'resp3'    % intrument response parameter 3 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(5,5) = tmp.val;
            
        case 'resp4'    % intrument response parameter 4 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(6,1) = tmp.val;
            
        case 'resp5'    % intrument response parameter 5 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(6,2) = tmp.val;
            
        case 'resp6'    % intrument response parameter 6 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(6,3) = tmp.val;
            
        case 'resp7'    % intrument response parameter 7 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(6,4) = tmp.val;
            
        case 'resp8'    % intrument response parameter 8 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(6,5) = tmp.val;
            
        case 'resp9'    % intrument response parameter 9 [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(7,1) = tmp.val;
            %}
            
        case 'stla'    % station latitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(7,2) = tmp.val;
            
        case 'stlo'    % station longitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(7,3) = tmp.val;
            
        case 'stel'    % station elevation
            assertions('number',tmp.hdr,tmp.val);
            h.float32(7,4) = tmp.val;
            
            %{
        case 'stdp'    % station depth [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(7,5) = tmp.val;
            %}
            
        case 'evla'    % event latitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(8,1) = tmp.val;
            
        case 'evlo'    % event longitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(8,2) = tmp.val;
            
            %{
        case 'evel'    % event elevation [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.float32(8,3) = tmp.val;
            %}
            
        case 'evdp'    % event depth
            assertions('number',tmp.hdr,tmp.val);
            h.float32(8,4) = tmp.val;
            
        case 'mag'    % event magnitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(8,5) = tmp.val;
            
        case 'user0'    % user-defined variable 0
            assertions('number',tmp.hdr,tmp.val);
            h.float32(9,1) = tmp.val;
            
        case 'user1'    % user-defined variable 1
            assertions('number',tmp.hdr,tmp.val);
            h.float32(9,2) = tmp.val;
            
        case 'user2'    % user-defined variable 2
            assertions('number',tmp.hdr,tmp.val);
            h.float32(9,3) = tmp.val;
            
        case 'user3'    % user-defined variable 3
            assertions('number',tmp.hdr,tmp.val);
            h.float32(9,4) = tmp.val;
            
        case 'user4'    % user-defined variable 4
            assertions('number',tmp.hdr,tmp.val);
            h.float32(9,5) = tmp.val;
            
        case 'user5'    % user-defined variable 5
            assertions('number',tmp.hdr,tmp.val);
            h.float32(10,1) = tmp.val;
            
        case 'user6'    % user-defined variable 6
            assertions('number',tmp.hdr,tmp.val);
            h.float32(10,2) = tmp.val;
            
        case 'user7'    % user-defined variable 7
            assertions('number',tmp.hdr,tmp.val);
            h.float32(10,3) = tmp.val;
            
        case 'user8'    % user-defined variable 8
            assertions('number',tmp.hdr,tmp.val);
            h.float32(10,4) = tmp.val;
            
        case 'user9'    % user-defined variable 9
            assertions('number',tmp.hdr,tmp.val);
            h.float32(10,5) = tmp.val;
            
        case 'dist'    % source receiver distance (km)
            assertions('number',tmp.hdr,tmp.val);
            h.float32(11,1) = tmp.val;
            
        case 'az'    % azimuth
            assertions('number',tmp.hdr,tmp.val);
            h.float32(11,2) = tmp.val;
            
        case 'baz'    % back azimuth
            assertions('number',tmp.hdr,tmp.val);
            h.float32(11,3) = tmp.val;
            
        case 'gcarc'    % great circle distance (deg)
            assertions('number',tmp.hdr,tmp.val);
            h.float32(11,4) = tmp.val;
            
            %{
        case 'internal2'    % SAC internal variable
            assertions('number',tmp.hdr,tmp.val);
            h.float32(11,5) = tmp.val;
            
        case 'internal3'    % SAC internal variable
            assertions('number',tmp.hdr,tmp.val);
            h.float32(12,1) = tmp.val;
            %}
            
        case 'depmen'    % mean amplitude
            assertions('number',tmp.hdr,tmp.val);
            h.float32(12,2) = tmp.val;
            
        case 'cmpaz'    % component azimuth
            if ischar(tmp.val), tmp.val = sac_ascii_azinc('cmpaz', tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.float32(12,3) = tmp.val;
            
        case 'cmpinc'    % component incident angle
            if ischar(tmp.val), tmp.val = sac_ascii_azinc('cmpinc', tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.float32(12,4) = tmp.val;
            
            %{
        case 'xminimum'    % minimum value of x (spectral files only)
            assertions('number',tmp.hdr,tmp.val);
            h.float32(12,5) = tmp.val;
            
        case 'xmaximum'    % maximum value of x (spectral files only)
            assertions('number',tmp.hdr,tmp.val);
            h.float32(13,1) = tmp.val;
            
        case 'yminimum'    % minimum value of y (spectral files only)
            assertions('number',tmp.hdr,tmp.val);
            h.float32(13,2) = tmp.val;
            
        case 'ymaximum'    % maximum value of y (spectral files only)
            assertions('number',tmp.hdr,tmp.val);
            h.float32(13,3) = tmp.val;
            
        case 'unused1'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(13,4) = tmp.val;
            
        case 'unused2'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(13,5) = tmp.val;
            
        case 'unused3'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(14,1) = tmp.val;
            
        case 'unused4'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(14,2) = tmp.val;
            
        case 'unused5'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(14,3) = tmp.val;
            
        case 'unused6'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(14,4) = tmp.val;
            
        case 'unused7'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.float32(14,5) = tmp.val;
            %}
            
        case 'nzyear'    % GMT year
            assertions('number',tmp.hdr,tmp.val);
            h.int32(1,1) = tmp.val;
            
        case 'nzjday'    % GMT julian day
            assertions('number',tmp.hdr,tmp.val);
            h.int32(1,2) = tmp.val;
            
        case 'nzhour'    % GMT hour
            assertions('number',tmp.hdr,tmp.val);
            h.int32(1,3) = tmp.val;
            
        case 'nzmin'    % GMT minute
            assertions('number',tmp.hdr,tmp.val);
            h.int32(1,4) = tmp.val;
            
        case 'nzsec'    % GMT second
            assertions('number',tmp.hdr,tmp.val);
            h.int32(1,5) = tmp.val;
            
        case 'nzmsec'    % GMT millisecond
            assertions('number',tmp.hdr,tmp.val);
            h.int32(2,1) = tmp.val;
            
        case 'nvhdr'    % header version number
            assertions('number',tmp.hdr,tmp.val);
            h.int32(2,2) = tmp.val;
            
        case 'norid'    % origin id
            assertions('number',tmp.hdr,tmp.val);
            h.int32(2,3) = tmp.val;
            
        case 'nevid'    % event id
            assertions('number',tmp.hdr,tmp.val);
            h.int32(2,4) = tmp.val;
            
        case 'npts'    % number of data points
            assertions('number',tmp.hdr,tmp.val);
            h.int32(2,5) = tmp.val;
            
            %{
        case 'internal4'    % SAC internal variable
            assertions('number',tmp.hdr,tmp.val);
            h.int32(3,1) = tmp.val;
            %}
            
        case 'nwfid'    % waveform id
            assertions('number',tmp.hdr,tmp.val);
            h.int32(3,2) = tmp.val;
            
            %{
        case 'nxsize'    % spectral length (spectral files only)
            assertions('number',tmp.hdr,tmp.val);
            h.int32(3,3) = tmp.val;
            
        case 'nysize'    % spectral width (spectral files only)
            assertions('number',tmp.hdr,tmp.val);
            h.int32(3,4) = tmp.val;
            
        case 'unused8'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(3,5) = tmp.val;
            %}
            
        case 'iftype'    % type of file
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(4,1) = tmp.val;
            
        case 'idep'    % type of independent variable
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(4,2) = tmp.val;
            
        case 'iztype'    % reference time equivalence
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(4,3) = tmp.val;
            
            %{
        case 'unused9'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(4,4) = tmp.val;
            
        case 'iinst'    % type of recording instrument [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.int32(4,5) = tmp.val;
            
        case 'istreg'    % station geographic region [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.int32(5,1) = tmp.val;
            
        case 'ievreg'    % event geographic region [not currently used]
            assertions('number',tmp.hdr,tmp.val);
            h.int32(5,2) = tmp.val;
            %}
            
        case 'ievtyp'    % type of event
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(5,3) = tmp.val;
            
            %{
        case 'iqual'    % quality of data [not currently used]
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(5,4) = tmp.val;
            
        case 'isynth'    % synthetic data flag [not currently used]
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(5,5) = tmp.val;
            %}
            
        case 'imagtyp'    % magnitude type
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(6,1) = tmp.val;
            
        case 'imagsrc'    % source of magnitude information
            % turn a string header value into an integer (e.g., 'ib' -> 9)
            if ischar(tmp.val), tmp.val = sac_ascii_value(tmp.val); end
            assertions('number',tmp.hdr,tmp.val);
            h.int32(6,2) = tmp.val;
            
            %{
        case 'unused10'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(6,3) = tmp.val;
            
        case 'unused11'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(6,4) = tmp.val;
            
        case 'unused12'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(6,5) = tmp.val;
            
        case 'unused13'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(7,1) = tmp.val;
            
        case 'unused14'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(7,2) = tmp.val;
            
        case 'unused15'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(7,3) = tmp.val;
            
        case 'unused16'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(7,4) = tmp.val;
            
        case 'unused17'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(7,5) = tmp.val;
            %}
            
        case 'leven'    % TRUE if data is evenly spaced
            assertions('logical',tmp.hdr,tmp.val);
            h.int32(8,1) = tmp.val;
            
        case 'lpspol'    % TRUE if station components have a positive polarity
            assertions('logical',tmp.hdr,tmp.val);
            h.int32(8,2) = tmp.val;
            
        case 'lovrok'    % TRUE if it is okay to overwrite this file on disk.
            assertions('logical',tmp.hdr,tmp.val);
            h.int32(8,3) = tmp.val;
            
        case 'lcalda'    % TRUE if DIST, AZ, BAZ, and GCARC are to be calculated from station and event coordinates.
            assertions('logical',tmp.hdr,tmp.val);
            h.int32(8,4) = tmp.val;
            
            %{
        case 'unused17'    % SAC unused header space
            assertions('number',tmp.hdr,tmp.val);
            h.int32(8,5) = tmp.val;
            %}
            
        case 'kstnm'    % station name
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(1,1:8) = tmp.val;
            
        case 'kevnm'    % event name
            assertions('string16',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-16s',tmp.val);  % pad end w/spaces
            h.char(1,9:24) = tmp.val;
            
        case 'ko'    % event origin time string
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(2,9:16) = tmp.val;
            
        case 'ka'    % first arrival time string
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(2,17:24) = tmp.val;
            
        case 'kt0'    % user-defined pick string 0
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(3,1:8) = tmp.val;
            
        case 'kt1'    % user-defined pick string 1
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(3,9:16) = tmp.val;
            
        case 'kt2'    % user-defined pick string 2
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(3,17:24) = tmp.val;
            
        case 'kt3'    % user-defined pick string 3
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(4,1:8) = tmp.val;
            
        case 'kt4'    % user-defined pick string 4
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(4,9:16) = tmp.val;
            
        case 'kt5'    % user-defined pick string 5
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(4,17:24) = tmp.val;
            
        case 'kt6'    % user-defined pick string 6
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(5,1:8) = tmp.val;
            
        case 'kt7'    % user-defined pick string 7
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(5,9:16) = tmp.val;
            
        case 'kt8'    % user-defined pick string 8
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(5,17:24) = tmp.val;
            
        case 'kt9'    % user-defined pick string 9
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(6,1:8) = tmp.val;
            
        case 'kf'    % end of event time string
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(6,9:16) = tmp.val;
            
        case 'kuser0'    % user-defined variable string 0
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(6,17:24) = tmp.val;
            
        case 'kuser1'    % user-defined variable string 1
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(7,1:8) = tmp.val;
            
        case 'kuser2'    % user-defined variable string 2
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(7,9:16) = tmp.val;
            
        case 'kcmpnm'    % component name
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(7,17:24) = tmp.val;
            
        case 'knetwk'    % name of seismic network
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(8,1:8) = tmp.val;
            
        case 'kinst'    % generic name of instrument
            assertions('string8',tmp.hdr,tmp.val);
            tmp.val = sprintf('%-8s',tmp.val);  % pad end w/spaces
            h.char(8,17:24) = tmp.val;
            
        otherwise
            error('ERROR:  "%s" is not a supported SAC header.', tmp.hdr);
            
    end
    
    
end

% open a new sac file for writing
fid = fopen(filename,'w');

% make sure we're at the beginning of the file
fseek(fid,0,'bof');

% write the sac file
fwrite(fid,h.float32','float32'); % float32 headers
fwrite(fid,h.int32','int32');     % int32 (and logical) headers
fwrite(fid,h.char','char');       % char headers
fwrite(fid,waveform,'float32');   % waveform data

% close and finish writing
fclose(fid);

% display success
fprintf('Wrote to file: %s\n', filename);


end


function val = sac_ascii_azinc(string, value)
% this function calculates the correct values for azimuth and incident
% angle given ascii inputs of 'Z', 'N', or 'E' instead of numbers

switch string
    case 'cmpaz'  % component azimuth
        switch lower( value )
            case 'z'  % if this is the z component
                val = 0;
            case 'n'  % if this is the n component
                val = 0;
            case 'e'  % if this is the e component
                val = 90;
            otherwise
                error(['ERROR:  "%s" not a supported component name. \n' ...
                    'Please try "Z", "N", "E", or directly input degrees for the "CMPAZ" header.'], value);
        end
    case 'cmpinc'  % component incident angle
        switch lower( value )
            case 'z'  % if this is the z component
                val = 0;
            case 'n'  % if this is the n component
                val = 90;
            case 'e'  % if this is the e component
                val = 90;
            otherwise
                error(['ERROR:  "%s" not a supported component name. \n' ...
                    'Please try "Z", "N", "E", or directly input degrees for the "CMPINC" header.'], value);
        end
end

end


function val = sac_ascii_value(string)
% this function takes strings as an argument and returns a value based on
% http://www.iris.edu/manuals/sac/SAC_Manuals/FileFormatPt2.html

switch lower( string )
    case 'itime',    val = 01;  % IFTYPE Time series file
    case 'irlim',    val = 02;  % IFTYPE Spectral file--real and imaginary
    case 'iamph',    val = 03;  % IFTYPE Spectral file--amplitude and phase
    case 'ixy',      val = 04;  % IFTYPE General x versus y data
    case 'iunkn',    val = 05;  % IDEP, IZTYPE, IEVTYP Unknown
    case 'idisp',    val = 06;  % IDEP Displacement in nm
    case 'ivel',     val = 07;  % IDEP Velocity in nm/sec
    case 'iacc',     val = 08;  % IDEP Velocity in nm/sec/sec
    case 'ib',       val = 09;  % IZTYPE Begin time
    case 'iday',     val = 10;  % IZTYPE Midnight of reference GMT day
    case 'io',       val = 11;  % IZTYPE Event origin time
    case 'ia',       val = 12;  % IZTYPE First arrival time
    case 'it0',      val = 13;  % IZTYPE User defined time pick 0
    case 'it1',      val = 14;  % IZTYPE User defined time pick 1
    case 'it2',      val = 15;  % IZTYPE User defined time pick 2
    case 'it3',      val = 16;  % IZTYPE User defined time pick 3
    case 'it4',      val = 17;  % IZTYPE User defined time pick 4
    case 'it5',      val = 18;  % IZTYPE User defined time pick 5
    case 'it6',      val = 19;  % IZTYPE User defined time pick 6
    case 'it7',      val = 20;  % IZTYPE User defined time pick 7
    case 'it8',      val = 21;  % IZTYPE User defined time pick 8
    case 'it9',      val = 22;  % IZTYPE User defined time pick 9
    case 'iradnv',   val = 23;  % undocumented
    case 'itannv',   val = 24;  % undocumented
    case 'iradev',   val = 25;  % undocumented
    case 'itanev',   val = 26;  % undocumented
    case 'inorth',   val = 27;  % undocumented
    case 'ieast',    val = 28;  % undocumented
    case 'ihorza',   val = 29;  % undocumented
    case 'idown',    val = 30;  % undocumented
    case 'iup',      val = 31;  % undocumented
    case 'illlbb',   val = 32;  % undocumented
    case 'iwwsn1',   val = 33;  % undocumented
    case 'iwwsn2',   val = 34;  % undocumented
    case 'ihglp',    val = 35;  % undocumented
    case 'isro',     val = 36;  % undocumented
    case 'inucl',    val = 37;  % IEVTYP Nuclear event
    case 'ipren',    val = 38;  % IEVTYP Nuclear pre-shot event
    case 'ipostn',   val = 39;  % IEVTYP Nuclear post-shot event
    case 'iquake',   val = 40;  % IEVTYP Earthquake
    case 'ipreq',    val = 41;  % IEVTYP Foreshock
    case 'ipostq',   val = 42;  % IEVTYP Aftershock
    case 'ichem',    val = 43;  % IEVTYP Chemical explosion
    case 'iother',   val = 44;  % IQUAL, IEVTYP Other
    case 'igood',    val = 45;  % IQUAL Good data
    case 'iglch',    val = 46;  % IQUAL Glitches
    case 'idrop',    val = 47;  % IQUAL Dropouts
    case 'ilowsn',   val = 48;  % IQUAL Low signal to noise ratio
    case 'irldta',   val = 49;  % ISYNTH Real data
    case 'ivolts',   val = 50;  % IDEP Velocity in volts
    case 'ixyz',     val = 51;  % IFTYPE General XYZ (3-D) file
    case 'imb',      val = 52;  % IMAGTYP Bodywave Magnitude
    case 'ims',      val = 53;  % IMAGTYP Surfacewave Magnitude
    case 'iml',      val = 54;  % IMAGTYP Local Magnitude
    case 'imw',      val = 55;  % IMAGTYP Moment Magnitude
    case 'imd',      val = 56;  % IMAGTYP Duration Magnitude
    case 'imx',      val = 57;  % IMAGTYP User Defined Magnitude
    case 'ineic',    val = 58;  % IMAGSRC National Earthquake Information Center
    case 'ipde',     val = 59;  % IMAGSRC Preliminary Determination of Epicenter
    case 'iisc',     val = 60;  % IMAGSRC Internation Seismological Centre
    case 'ireb',     val = 61;  % IMAGSRC Reviewed Event Bulletin
    case 'iusgs',    val = 62;  % IMAGSRC US Geological Survey
    case 'ibrk',     val = 63;  % IMAGSRC UC Berkeley
    case 'icaltech', val = 64;  % IMAGSRC California Institute of Technology
    case 'illnl',    val = 65;  % IMAGSRC Lawrence Livermore National Laboratory
    case 'ievloc',   val = 66;  % IMAGSRC Event Location (computer program)
    case 'ijsop',    val = 67;  % IMAGSRC Joint Seismic Observation Program
    case 'iuser',    val = 68;  % IMAGSRC The individual using SAC2000
    case 'iunknown', val = 69;  % IMAGSRC unknown
    case 'iqb',      val = 70;  % IEVTYP Quarry or mine blast confirmed by quarry
    case 'iqb1',     val = 71;  % IEVTYP Quarry/mine blast with designed shot info-ripple fired
    case 'iqb2',     val = 72;  % IEVTYP Quarry/mine blast with observed shot info-ripple fired
    case 'iqbx',     val = 73;  % IEVTYP Quarry or mine blast - single shot
    case 'iqmt',     val = 74;  % IEVTYP Quarry/mining-induced events: tremors and rockbursts
    case 'ieq',      val = 75;  % IEVTYP Earthquake
    case 'ieq1',     val = 76;  % IEVTYP Earthquakes in a swarm or aftershock sequence
    case 'ieq2',     val = 77;  % IEVTYP Felt earthquake
    case 'ime',      val = 78;  % IEVTYP Marine explosion
    case 'iex',      val = 79;  % IEVTYP Other explosion
    case 'inu',      val = 80;  % IEVTYP Nuclear explosion
    case 'inc',      val = 81;  % IEVTYP Nuclear cavity collapse
    case 'io_',      val = 82;  % IEVTYP Other source of known origin
    case 'il',       val = 83;  % IEVTYP Local event of unknown origin
    case 'ir',       val = 84;  % IEVTYP Regional event of unknown origin
    case 'it',       val = 85;  % IEVTYP Teleseismic event of unknown origin
    case 'iu',       val = 86;  % IEVTYP Undetermined or conflicting information
    otherwise
        error(['ERROR: "%s" is not a known header value.\n' ...
            '  Check the SAC manual for proper header values.\n' ...
            '  http://www.iris.edu/manuals/sac/SAC_Manuals/FileFormatPt2.html'],string);
end

return

end

function assertions(type,hdr,val)
% This function saves space so that repeated assertions don't have to be
% copied and pasted a thousand times in the code.

switch type
    case 'float32'
        assert( isequal( size(val), [14 5] ), ...  % make sure that the header block is the correct size )
            ['ERROR: "%s" is not the proper size of the 14x5float32 header block.\n' ...
            'Please see http://www.iris.edu/software/sac/manual/file_format.html'], hdr);
        
    case 'int32'
        assert( isequal( size(val), [8 5] ), ...  % make sure that the header block is the correct size )
            ['ERROR: "%s" is not the proper size of the 8x5 int32 header block.\n' ...
            'Please see http://www.iris.edu/software/sac/manual/file_format.html'], hdr);
        
    case 'char'
        assert( isequal( size(val), [8 24] ), ...  % make sure that the header block is the correct size )
            ['ERROR: "%s" is not the proper size of the 8x24 char header block.\n' ...
            'Please see http://www.iris.edu/software/sac/manual/file_format.html'], hdr);
        
    case 'number'
        assert( isnumeric(val) & isscalar(val), ...  % make sure argument is a number
            'ERROR: "%s" should be a 1x1 real matrix.', hdr);
        
    case 'logical'
        assert( islogical(val) | isequal(val,1) | isequal(val,0) & isscalar(val), ...  % and logical or 0 or 1
            'ERROR: "%s" should be a 1x1 logical matrix (i.e., 0 or 1).', hdr);
        
    case 'string8'
        assert( ischar(val), 'ERROR: "%s" must be a string.', hdr);  % make sure the argument is a string
        assert( numel(val) <= 8, 'ERROR: "%s" must be <= 8 characters long', hdr);  % argument length <= 8 characters
        
    case 'string16'
        assert( ischar(val), 'ERROR: "%s" must be a string.', hdr);  % make sure the argument is a string
        assert( numel(val) <= 16, 'ERROR: "%s" must be <= 16 characters long', hdr);  % argument length <= 16 characters
end

end

