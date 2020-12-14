function seminar2(subj, rec) 
  db = "baza";
  subject = string(subj);
  record=[];
  for i=0:2
      record = [record string(num2str(rec+(i*4), '%02d'))];
  end
  t1s = {};
  t2s = {};

  for i=1:size(record,2)
    if (strcmp(db,"")==0) 
      recName=strcat("/",db,"/",subject,"/",subject,"R",record(:,i),".edf");
    else
      recName=strcat(subject,'R',record(i),'.edf');
    end
    recName=convertStringsToChars(recName);
    disp(recName);
    [sig, fs, tm] = rdsamp(recName, 1:64); % 
    [T0, T1, T2] = getIntervals(recName,'event', fs, size(sig,1));
    sig=sig';

    sec = 4.0;
    fd = fs* sec;
    
    [m,n] = size(sig);
    
    for j=1:size(T1, 1)
      if (T1(j,1) + fd) > n
          t1s{end+1}=sig(:, T1(j,1):n);
      else
          t1s{end+1}=sig(:, T1(j,1):T1(j,1) + fd);
      end
    end
disp('T1')
disp(T1) % Tu so zacetni in koncni indeksi intervalov zamisljanja leve roke

    for j=1:size(T2,1)
      if (T2(j,1) + fd) > n
          t2s{end+1}=sig(:, T2(j,1):n);
      else
          t2s{end+1}=sig(:, T2(j,1):T2(j,1) + fd);
      end
      
    end
    
    disp('T2')
disp(T2)
%disp(T2) % Tu so zacetni in koncni indeksi intervalov zamisljanja desne roke
  end
  
  
if size(t1s,2) > 0 && size(t2s,2) > 0 %varovalka, če subjekt nima intervalov zamišljanja T1 ali T2
    
    
N_int_povp = 5; % Stevilo intervalov za povprecna intervala (pet)
% Povprecje za intervale leve roke
t1s_mean = cell2mat(t1s(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2: N_int_povp  % size(T1,2)  % V povprecje vkljucim prvih pet intervalov leve roke, size(T1,2) bi bili vsi intervali
  t1s_mean = t1s_mean + cell2mat(t1s(j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1s_mean = t1s_mean / size(T1,2);
size(t1s_mean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Povprecje za intervale desne roke
t2s_mean = cell2mat(t2s(1));
for j = 2: N_int_povp  % size(T2,2) % V povprecje vkljucim prvih pet intervalov desne roke, size(Ts,2) bi bili vsi intervali desno
  t2s_mean = t2s_mean + cell2mat(t2s(j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2s_mean = t2s_mean / size(T2,2);
size(t2s_mean)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp1 = t1s_mean;
tmp2 = t2s_mean;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Narisemo povprecna intervala zamisljanja leve in desne roke
figure(603)
for t1=1:2
%plot(tmp1(1,:))
    plot(tmp1(t1,:)) % Prvi signal, Slika je v samplih signala po X osi
    ylim([-1000, 1000]);
    hold on                        % naslednji signa bo v isti sliki
    %plot(tmp1(size(tmp1,1),:))       % Drugi signal, Slika je v samplih signala po X osi
%ylim([-1000, 1000]);
end

figure(604)
for t2=1:2
%plot(tmp2(1,:))
    plot(tmp2(t2,:)) % Prvi signal, Slika je v samplih signala po X osi
    ylim([-1000, 1000]);
    hold on                        % naslednji signa bo v isti sliki
    %plot(tmp2(size(tmp2,1),:))       % Drugi signal, Slika je v samplih signala po X osi
%ylim([-1000, 1000]);
end
%% fs = 160

%%%%  f = [0 8 8 13 13 fs/2]/(fs/2); 

% Bandpass filter from 8 to 13 Hz
  f = [0  7  8  13  14  fs/2]/(fs/2); 
  a = [0  0  1   1   0  0 ];
  n= 35; 
  b = firls (n, f, a);

%  Filtriranje je vkluceno kasneje
%   [W] = f_CSP(filter(b, 1, cell2mat(t1s(1))), filter(b, 1, cell2mat(t2s(1))));

% Metoda CSP Common Spatial Patterns
%   [W] = f_CSP( cell2mat(t1s(1)),  cell2mat(t2s(1)));
% Metoda CSP Common Spatial Patterns

[W] = f_CSP( t1s_mean,  t2s_mean);


  %t1s(1)=[];
  %t2s(1)=[];
  
  lvt1=[];
  lvt2=[];

%%%%

lengths_of_intervals_in_samlpes_left = [];
lengths_of_intervals_in_seconds_left = [];
%%%%

  for i=1:size(t1s,2)
      % tmp = (W*filter(b, 1,cell2mat(t1s(i))));

% Signali se prevedejo v prostor stanj 64 -> 64,
% zanimata nas le prvi (postane PRVI signal v prostoru stanj) 
% in 64-ti (postane DRUGI signal v prostoru stanj)

        tmp = (W*cell2mat(t1s(i)));

% disp(length(tmp(1,:)));
% disp(length(tmp(64,:)));
% pause

% Slike signalov za intervale zamisljanja leve roke pred filtriranjem (slike od 1 naprej)

%%%%    tmp(size(tmp,1),:)       drugi signal v prostoru stanj (na polozaju 64 v matriki) pred filtriranjem

lengths_of_intervals_in_samlpes_left = [lengths_of_intervals_in_samlpes_left, length(tmp(1,:))];
lengths_of_intervals_in_seconds_left = [lengths_of_intervals_in_seconds_left, length(tmp(1,:))*1/fs];

% Tu se filtrirata prvi in drugi signal v prostoru stanj
      tmp = [tmp(1,:).', tmp(size(tmp,1),:).' ].';
      tmp = filter(b, 1, tmp);



%%%%     tmp(2,:)                drugi signal v prostoru stanj PO filtriranju (na polozaju 2 v matriki)

%% Na tem mestu se filtrirani intervali zamisljanj leve roke (po dva signala v prostoru stanj)
%% lahko zapisejo na disk

% Izracun znacilk za intervale zamisljanja aktivnosti leve roke
      lvt1(i,1) = log(var(tmp(1,:)));
      lvt1(i,2) = log(var(tmp(2,:)));     
  end

% Dolzine intervalov zamisljanja aktivnosti leve roke v vzorcih
lengths_of_intervals_in_samlpes_left
% Dolzine intervalov zamisljanja aktivnosti leve roke v sekundah
lengths_of_intervals_in_seconds_left

%%%%

lengths_of_intervals_in_samlpes_right = [];
lengths_of_intervals_in_seconds_right = [];

  for i=1:size(t2s,2)
      % tmp = (W*filter(b, 1,cell2mat(t2s(i))));

% Signali se prevedejo v prostor stanj 64 -> 64,
% zanimata nas le prvi (postane PRVI signal v prostoru stanj) 
% in 64-ti (postane DRUGI signal v prostoru stanj)

        tmp = (W*cell2mat(t2s(i)));

% disp(length(tmp(1,:)));
% disp(length(tmp(64,:)));
% pause

%%%%    tmp(size(tmp,1),:)       drugi signal v prostoru stanj (na polozaju 64 v matriki) pred filtriranjem

lengths_of_intervals_in_samlpes_right = [lengths_of_intervals_in_samlpes_right, length(tmp(1,:))];
lengths_of_intervals_in_seconds_right = [lengths_of_intervals_in_seconds_right, length(tmp(1,:))*1/fs];

% Tu se filtrirata prvi in drugi signal v prostoru stanj
      tmp = [tmp(1,:).', tmp(size(tmp,1),:).' ].';
      tmp = filter(b, 1, tmp);

%%%%    tmp(2,:)                 drugi signal v prostoru stanj PO filtriranju (na polozaju 2 v matriki)


% Izracun znacilk za intervale zamisljanja aktivnosti desne roke roke
      lvt2(i,1) = log(var(tmp(1,:)));
      lvt2(i,2) = log(var(tmp(2,:)));     
  end

% Dolzine intervalov zamisljanja aktivnosti desne roke v vzorcih
lengths_of_intervals_in_samlpes_right
% Dolzine intervalov zamisljanja aktivnosti desne roke v sekundah
lengths_of_intervals_in_seconds_right

%%%



figure(500) 

  scatter(lvt1(:,1), lvt1(:,2));
  hold on
  scatter(lvt2(:,1), lvt2(:,2)); 
  
      featVFile = strcat(subject,'featureVectors.txt');
      classFile = strcat(subject,'referenceClass.txt');

      fvf = fopen(featVFile, "wt");
      rcf = fopen(classFile, "wt");
  
  for i=1:size(lvt1,1)
      fprintf(fvf, "%.8f %.8f\n", lvt1(i,1), lvt1(i,2));
      fprintf(rcf, "T1\n");
  end
  
  for i=1:size(lvt2,1)
      fprintf(fvf, "%.8f %.8f\n", lvt2(i,1), lvt2(i,2));
      fprintf(rcf, "T2\n");
  end
  fclose(fvf);
  fclose(rcf);
      
end

end

function [newVectors, meanValue] = remmean(vectors)
  newVectors = zeros (size (vectors));
  meanValue = mean (vectors')';
  newVectors = vectors - meanValue * ones (1,size (vectors, 2));
end

function [T0,T1, T2] = getIntervals(recName, anotator, fs, siglen)
% function returns arrays of interval start and stop times in samples 
% for the three different interval types; the meaning of annotations (T1,
% T2, T3) is determined according to the record number;
% input arguments:
%   - recName: string containing record name e.g 'S001R03.edf'
%   - anotator: string containing annotator name e.g. 'event'
%   - fs: sampling frequency of the record e.g. 160
%   - siglen: length of a data file in samples
% output arguments:
%   - T0: array of start and end times in samples for interval T0
%   - T1: array of start and end times in samples for interval T1
%   - T2: array of start and end times in samples for interval T2
  T0=[];
  T1=[];
  T2=[];
  [a, t, s, c, n, cmt] = rdann(recName, anotator);
  cmt=string(cmt);
  for i=1:size(a,1)
    splt=strsplit(cmt(i), " ");
    dur=str2double(splt(3))*fs;
    if (i<size(a,1)) 
        stop=a(i+1)-1;
    else
        stop=a(i)+dur;
    end
    stop = min(stop, siglen);
    if (splt(1)=="T0") 
      T0=[T0; [a(i) stop]];
    elseif (splt(1)=="T1")
      T1=[T1; [a(i) stop]];
    elseif (splt(1)=="T2")
      T2=[T2; [a(i) stop]];
    end
  end
end

function [result] = f_CSP(varargin)
% This code is for calulating the projection matrix for CSP 
% Haider Raza, Intelligent System Research Center, University of Ulster, Northern Ireland, UK.
%     Raza-H@email.ulster.ac.uk
%     Date: 03-Oct-2014
% Input:
%             
%       left:  left hand data
%       right: right hand data
% 
% Output:
%        left: Left hand data
%        right: right hand data    

    if (nargin ~= 2)
        disp('Must have 2 classes for CSP!')
    end
    
    Rsum=0;
    %finding the covariance of each class and composite covariance
    for i = 1:nargin 
        %mean here?
        R{i} = ((varargin{i}*varargin{i}')/trace(varargin{i}*varargin{i}'));%instantiate me before the loop!
        %Ramoser equation (2)
        Rsum=Rsum+R{i};
    end
   
    %   Find Eigenvalues and Eigenvectors of RC
    %   Sort eigenvalues in descending order
    [EVecsum,EValsum] = eig(Rsum);
    [EValsum,ind] = sort(diag(EValsum),'descend');
    EVecsum = EVecsum(:,ind);
    
    %   Find Whitening Transformation Matrix - Ramoser Equation (3)
        W = sqrt(inv(diag(EValsum))) * EVecsum';
    
    
    for k = 1:nargin
        S{k} = W * R{k} * W'; %       Whiten Data Using Whiting Transform - Ramoser Equation (4)
    end
   
    %generalized eigenvectors/values
    [B,D] = eig(S{1},S{2});
    % Simultanous diagonalization
			% Should be equivalent to [B,D]=eig(S{1});
    
    [D,ind]=sort(diag(D));B=B(:,ind);
    
    %Resulting Projection Matrix-these are the spatial filter coefficients
    result = B'*W;

end

