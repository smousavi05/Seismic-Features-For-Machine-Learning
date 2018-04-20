% This function is made to compute the polarity parameters of 3 components
% seismogram. The adavnatge of this one is that it only computes polarity around
% user specified value since this analysis is really time consuming
% Input:    array_ind   array of indices around whom we want to compute
%                       polarity
%           window      size of window around indices
%           samplrate   sampling rate
%           trace?      traces
% Output:   rect
%           azi
%           dip

function [rectilinP,azimuthP,dipP]=fast_polar_analysis(array_ind,window,samplrate,tracex,tracey,tracez)



dt=1/samplrate;

ntr=1;

%%%% CAREFULLY
%%%% assign half-length of the moving time window
halfInterval=2; % in seconds
w= round(halfInterval/dt);


N_win=round(window*samplrate);

n_max=max([length(tracex) length(tracey) length(tracez)]);
ns=min([length(tracex) length(tracey) length(tracez)]);

ind_vec=[];
for i=1:length(array_ind)
    ind_vec=[ind_vec;(array_ind(i)-N_win+1:array_ind(i)+N_win)'];
end
ind_vec(ind_vec<=0)=[];
ind_vec(ind_vec>ns)=[];
ind_vec=sort(ind_vec);
ind_vec(diff(ind_vec)==0)=[];

% We have to be sure that the samples are ok with polarization window
% i.e first sample should be a least bigger than w so than that
% first_ind-w>=1

ind_vec(ind_vec-w<=0)=[];
ind_vec(ind_vec+w>ns)=[];


%%% which contain ?dt? the sampling time in seconds,
%%% ?ns? number of samples,
%%% ?ntr? number of traces,
%%% ?rcvx?,?rcvy?,?rcvz? receiver coordinates,
%%% ?tracex?,?tracey?,?tracez? the data

rectilinP=zeros(ns,ntr);
azimuthP=zeros(ns,ntr);
dipP=zeros(ns,ntr);

%%%% define start and end time for the analysis
%trstart=w+1;
%trend=ns-w;
%%% loop over traces
for h=1:ntr
%%%% loop over moving time window
    for k=[ind_vec']
        %clear xp yp zp MP pP DP
        xP=tracex(k-w:k+w,h);
        yP=tracey(k-w:k+w,h);
        zP=tracez(k-w:k+w,h);
        MP=cov([xP yP zP]);%covariance(xP,yP,zP);
        [pP,DP] = eig(MP);
        %%% DP contains the eigenvalues of the covariance matrix, with
        %%% DP(3,3)>DP(2,2)>DP(1,1)
        %%% pP contains the eigenvectors, where the first column is the
        %%% eigenvector that corresponds to the smallest eigenvalue, the
        %%% second one to the intermedian eigenvalue and the third one to
        %%% the largest eigenvalue (this one shows the dominant particle motion)
        rectilinP(k,h)=1-((DP(1,1)+DP(2,2))/(2*DP(3,3)));
        azimuthP(k,h)=atan(pP(2,3)/pP(1,3))*180/pi;
        dipP(k,h)=atan(pP(3,3)/sqrt(pP(2,3)^2+pP(1,3)^2))*180/pi;
    end;
end;

        % Shift those value by window size

%         rectilinP=[NaN(w,h); rectilinP];
%         azimuthP=[NaN(w,h); azimuthP];
%         dipP=[NaN(w,h); dipP];
end