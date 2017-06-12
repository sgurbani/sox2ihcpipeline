%%Tissue Area Count script
%%Date: June 8th, 2015
%
%%This script will read in an image file specified by IMG_FILE
%
% Input params:
%   img_name - name of the image file (no path, no extension)
%   input_dir - location of the image
%   output_dir - location to save area calculations (as img_name_area.txt)
%                   and TIFFs (as img_name_Level_Resolution.jpg)

function TissueAreaCount( img_name , input_dir, output_dir)
tic;

%% START PARAMS

img_file = [input_dir filesep img_name '.ndpi'];

disp(img_file);

%histogram params
bins = 1:255;

%Gaussian params
wLarge = 8; %large Gaussian window for filtering
wSmall = 4; %small Guassian window for filtering

%glass estimate
MuGlass = 20;
dMuGlass = 20;
SigSqGlass = 4;

%tissue mean estimate
MuTissue = 110;
dMuTissue = 30;
SigSqTissue = 10;

%nuclei mean estimate
MuNuclei = 200;
dMuNuclei = 40;
SigSqNuclei = 15;

%% END PARAMS

fprintf('Starting image %s.\n', img_name');

%% BUILD MODEL

fprintf('Building mixture model...');

%add path to subroutines
addpath('/home/lcoop22/ImageAnalysis/IHCPipeline/');         
addpath('/home/lcoop22/ImageAnalysis/OpenSlide/');                  %Openslide
addpath('/home/sgurban/ImageAnalysis/ColorDeconvolution/');         %Color deconvolution
addpath('/home/sgurban/Bioinformatics/MixtureModeling/');           %Mixture Modeling
addpath('/home/lcoop22/ImageAnalysis/MaxFlowBinarization/');        %graph cuts
load('/home/sgurban/ImageAnalysis/ColorDeconvolution/H.DAB.mat');   %stain matrices

%see if we can open the image file
if ~openslide_can_open( img_file )
    disp('Error - could not open image');
    return;
end

%get level params and dimensions
[dims, factors, objective, MppX, MppY] = openslide_check_levels( img_file );

%find the factor with approximately 1k height (usually about 2k width)
ind_of_1k_res = find( dims(:,2) < 2000, 1);

%open this file
img_1k = openslide_read_regions( img_file, ind_of_1k_res - 1, 0, 0, dims(ind_of_1k_res,1), dims(ind_of_1k_res,2) );
img_1k = img_1k{1}; %get only element of cell array
img_1k = permute(img_1k, [2 1 3]);   %transpose for viewability

%use mode to get estimate of background level
%background_est = squeeze(mode(mode(img_1k)));

%remove background blocks from edge
background_inds = ( img_1k(:,:,1) == 0 ) & ( img_1k(:,:,2) == 0 ) & ( img_1k(:,:,3) == 0 );
edge_inds = imclearborder(background_inds);
background = xor(edge_inds, background_inds);

%separate hematoxylin and DAB channels
intensity = ColorDeconvolution( double(img_1k), M, [1 1 0], false);
scaled = uint8(intensity);
scaled = imfilter(imcomplement(scaled(:,:,1)), fspecial('gaussian', wLarge, (wLarge/2)/3));

%mask for tissue region
mask = scaled >= MuGlass + 2*SigSqGlass;
mask = imfill(mask, 'holes');
mask = imdilate(mask, strel('disk', 10));

%apply background border mask to this
mask = ~background & mask;

%make mask RGB
mask_rgb = uint8(zeros(size(img_1k)));
mask_rgb(:,:,1) = uint8(mask);
mask_rgb(:,:,2) = uint8(mask);
mask_rgb(:,:,3) = uint8(mask);

%generate histogram for tissue region
h_tissue = histc(scaled(mask).', bins);
h_tissue = h_tissue / sum(h_tissue);

%find peaks in the regions of the histogram for Tissue, Nuclei, and Glass
tissuePeak = max(h_tissue( (bins >= MuTissue - dMuTissue) & (bins <= MuTissue + dMuTissue) ));
nucleiPeak = max(h_tissue( (bins >= MuNuclei - dMuNuclei) & (bins <= MuNuclei + dMuNuclei) ));
glassPeak = max(h_tissue(bins < MuTissue - dMuTissue));

%build params for curvefit model
x0 = zeros(1,8);            %initial estimate vector, initialization
x0(1) = glassPeak / (glassPeak + tissuePeak + nucleiPeak);      %p(glass)
x0(2) = tissuePeak / (glassPeak + tissuePeak + nucleiPeak);     %p(tissue)
x0(3) = MuGlass;
x0(4) = SigSqGlass;
x0(5) = MuTissue;
x0(6) = SigSqTissue;
x0(7) = MuNuclei;
x0(8) = SigSqNuclei;

%build function for curvefit model
fun = @(x,xdata)x(1)*normpdf(xdata, x(3), x(4)) + ...
    x(2)*normpdf(xdata, x(5), x(6)) + ...
    (1 - x(1) - x(2))*normpdf(xdata, x(7), x(8));

%get model estimate
globalModel = lsqcurvefit(fun, x0, bins, h_tissue, ...
    [0 0 MuGlass-dMuGlass 0  MuTissue-dMuTissue 0  MuNuclei-dMuNuclei 0 ], ...
    [1 1 MuGlass+dMuGlass 10 MuNuclei+dMuNuclei 50 MuNuclei+dMuNuclei 50], ...
    optimset('Display','off'));

%get params for tissue
gProbTissue = globalModel(2);
gMuTissue = globalModel(5);
gSigSqTissue = globalModel(6);

%get params for glass
gProbGlass = globalModel(1);
gMuGlass = globalModel(3);
gSigSqGlass = globalModel(4);

%get params for nuclei
gMuNuclei = globalModel(7);
gSigSqNuclei = globalModel(8);

%define transition between glass and tissue
MLglass = max( bins(gProbGlass*normpdf(bins, gMuGlass, gSigSqGlass) > ...
    gProbTissue*normpdf(bins, gMuTissue, gSigSqTissue)) );

%clear variables that we don't need
clearvars background background_inds edge_inds img_1k intensity mask ...
    mask_rgb scaled tissue;

fprintf('done.\n');
%% GET AREA CALCULATIONS ON HIGH RES IMAGES

%create output file
outfile = fopen( [output_dir filesep img_name '_area.txt'], 'w');
fprintf(outfile, '%-5s\t%-20s\t%-20s\n', 'Level', 'Resolution', 'Area (mm^2)');

for i = ind_of_1k_res:-1:1
    %check what the dims are, and if too big, continue
    if 6*prod(dims(i,:)) >= 20*intmax
        fprintf('Image too big at level %u.\n', factors(i));
        break;
    end

    fprintf('Calculating tissue area at level %u...', factors(i));
    %now, select a higher resolution image and estimate tissue area
    ind_hr = i;
    img_hr = openslide_read_regions( img_file, ind_hr - 1, 0, 0, dims(ind_hr,1), dims(ind_hr,2) );
    %img_hr = openslide_read_regions( img_file, ind_hr - 1, 8192, 36864, 4096, 4096 );
    img_hr = img_hr{1};

    %subtract blocks of background on boundary
    background_inds = ( img_hr(:,:,1)==0 ) & ( img_hr(:,:,2)==0 ) & ( img_hr(:,:,3)==0 );
    edge_inds = imclearborder(background_inds);
    background = xor(edge_inds, background_inds); clearvars background_inds edge_inds;

    %separate channels
    intensity = ColorDeconvolution( double(img_hr), M, [1 1 0], false );
    scaled = uint8(intensity); clearvars intensity;
    smoothed = imfilter( imcomplement( scaled(:,:,1) ), ...
        fspecial('gaussian', wLarge, (wLarge/2)/3)); clearvars scaled;

    %make mask
    mask = smoothed > MLglass;
    mask( mask & background ) = false;

    %apply MaxFlowBinarization to separate tissue from glass using global model
    pGlass = gProbGlass * normpdf( bins, gMuGlass, gSigSqGlass );
    pGlass = pGlass / sum(pGlass);  %normalize
    pTissue = gProbTissue * normpdf( bins, gMuTissue, gSigSqTissue );
    pTissue = pTissue / sum(pTissue); %normalize

    %format of params is: Image, Background, Foreground, sigma, weight
    tissue = MaxFlowBinarization( smoothed, pGlass, pTissue, 5, 4);
    tissue = tissue == 1 & ~background;

    %get tissue area estimate in mm^2
    tissue_area_pixels = length(find(tissue));
    tissue_area_microns2 = tissue_area_pixels * factors(ind_hr)*factors(ind_hr) * str2double(MppX) * str2double(MppY);
    tissue_area_mm2 = tissue_area_microns2 / 1e6;

    %create output image (full color and mask)
    img_out = zeros(size(img_hr,1), 2*size(img_hr,2), 3, 'uint8');
    img_out(:,1:size(img_hr,2),:) = uint8(img_hr); %store full color image on left
    img_out(:,size(img_hr,2)+1:end,1) = 255*uint8(mask);   %mask in right panel and red channel
    imwrite(img_out, [output_dir filesep img_name '_' num2str(factors(i)) '_' num2str(dims(i,1)) 'x' num2str(dims(i,2)) '.jpg']);

    %print to outputfile
    fprintf(outfile, '%-5s\t%-20s\t%7.4f\n', num2str(factors(i)), [num2str(dims(i,1)) ' x ' num2str(dims(i,2))], tissue_area_mm2);

    %clear variables
    clearvars img_hr img_out background smoothed mask;

    fprintf('done.\n');
end

fclose(outfile);

fprintf('Tissue area calculation complete. Output saved to %s\n', output_dir);
toc;

end
