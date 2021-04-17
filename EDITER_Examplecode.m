%%%%
% Created on 16 April 2021 
%
%
% Author : Sai Abitha Srinivas
%
% EDITER ALGORITHM 
%
% Input : kspace data with 5 emi detectors. 
% Modify to change # of emi detectors with "Nc" and commenting out extra
% "dftmp" and "df#p" 
% Modify from 2D to 3D data by adding extra for loop for different slices
%%%%

close all;
clear all;
clc;

load('data_BBEMI_2D_brainslice') %2D data from MGH Portable MRI scanner with broadband EMI

%% image size
ncol = size(datafft, 1);
nlin = size(datafft, 2);
Nc = 5;
%%%%%%%%%%%%%%%%%%%%%%
%% Initial pass using single PE line ( Nw =1 )
ksz_col = 0; % Deltakx = 1
ksz_lin = 0; % Deltaky = 1

%% kernels across pe lines
kern_pe = zeros( Nc * (2*ksz_col+1) * (2*ksz_lin+1), nlin);
    
for clin = 1:nlin
    
    noise_mat = [];
    
    pe_rng = clin;
    
    df1p = padarray(datanoise_fft_1(:,pe_rng), [ksz_col ksz_lin]);
    df2p = padarray(datanoise_fft_2(:,pe_rng), [ksz_col ksz_lin]);
    df3p = padarray(datanoise_fft_3(:,pe_rng), [ksz_col ksz_lin]);
    df4p = padarray(datanoise_fft_4(:,pe_rng), [ksz_col ksz_lin]);
    df5p = padarray(datanoise_fft_5(:,pe_rng), [ksz_col ksz_lin]);
        
    for col_shift = [-ksz_col:ksz_col]
        for lin_shift = [-ksz_lin:ksz_lin]
            
            dftmp = circshift(df1p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df2p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df3p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df4p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
             dftmp = circshift(df5p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
        end
    end
      
    %%%%%%%%%%%
    %% kernels
    gmat = reshape(noise_mat, size(noise_mat,1) * size(noise_mat,2), size(noise_mat,3));
    
    init_mat_sub = datafft(:,pe_rng);
    kern = gmat \ init_mat_sub(:);
    
    kern_pe(:,clin) = kern;
        
end

%%%%%%%%%%
%% look at correlation between the kernels

kern_pe_normalized = zeros(size(kern_pe));
for clin = 1:nlin

    kern_pe_normalized(:,clin) = kern_pe(:,clin) / norm(kern_pe(:,clin));

end

kcor = kern_pe_normalized'  * kern_pe_normalized;

%% threshold 
kcor_thresh = abs(kcor) > 5e-1;

%% start with full set of lines
aval_lins = [1:nlin];

%% window stack
win_stack = cell(nlin, 1);

cwin = 1;
while (~isempty(aval_lins))

        clin = min(aval_lins);
        pe_rng = [clin:clin+max(find(kcor_thresh(clin, clin:end)))-1];
        
        win_stack{cwin} = pe_rng;
        aval_lins = sort(setdiff(aval_lins, pe_rng), 'ascend');
    
        cwin = cwin + 1;
end

%% drop the empty entries
win_stack = win_stack(1:cwin-1);
    
ksz_col = 7; % Deltakx
ksz_lin = 0; % Deltaky

%% solution kspace
gksp = zeros(ncol, nlin);

for cwin = 1:length(win_stack)
    
    noise_mat = [];
    
    pe_rng = win_stack{cwin};
    
    df1p = padarray(datanoise_fft_1(:,pe_rng), [ksz_col ksz_lin]);
    df2p = padarray(datanoise_fft_2(:,pe_rng), [ksz_col ksz_lin]);
    df3p = padarray(datanoise_fft_3(:,pe_rng), [ksz_col ksz_lin]);
    df4p = padarray(datanoise_fft_4(:,pe_rng), [ksz_col ksz_lin]);
    df5p = padarray(datanoise_fft_5(:,pe_rng), [ksz_col ksz_lin]);
    
    for col_shift = [-ksz_col:ksz_col]
        for lin_shift = [-ksz_lin:ksz_lin]
            
            dftmp = circshift(df1p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df2p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df3p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df4p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
            dftmp = circshift(df5p(:,:), [col_shift, lin_shift]);
            noise_mat= cat(3, noise_mat, dftmp(ksz_col+1:end-ksz_col, ksz_lin+1:end-ksz_lin));
            
        end
    end
    
    
    %%%%%%%%%%%
    %% kernels
    gmat = reshape(noise_mat, size(noise_mat,1) * size(noise_mat,2), size(noise_mat,3));
    
    init_mat_sub = datafft(:,pe_rng);
    kern = gmat \ init_mat_sub(:);
        
    %% put the solution back
    tosub = reshape(gmat * kern, ncol, length(pe_rng));
    gksp(:,pe_rng) = init_mat_sub - tosub;
    
end

corr_img_opt_toep = fftshift(fftn(fftshift(gksp)));

x_range =[150:350];
y_range = [1:101];

figure(100);
subplot(211);
uncorr = fftshift(fftn(fftshift(datafft)));
imagesc(flipud(rot90(abs(uncorr(x_range,y_range))))); axis equal; colormap gray;
axis tight;
title([' primary uncorrected']);

subplot(212);
imagesc(flipud(rot90(abs(corr_img_opt_toep(x_range,y_range))))); axis equal; colormap gray;
axis tight;  
title([' corrected with EDITER' ]);
