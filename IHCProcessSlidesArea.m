function IHCProcessSlides(Slide, Desired, T, Folder)
%Uses openslide library to process tiles in slide for IHC analysis.
%Paths must be set to ColorDeconvolution, NewSegmentation, PoissonBinarization,
%MaxFlowBinarization and BoundaryValidator.
%inputs:
%Slide - string, path and filename of input file.
%Desired - scalar, magnification desired for analysis (typically 20).
%T - scalar, tilesize for analysis resolution (typically 4096).
%Folder - string, path for outputs.

%parameters
ModelingMag = 5; %magnification used to order tiles in terms of entropy content
Mfile = '/home/lcoop22/ImageAnalysis/ColorDeconvolution/H.DAB.mat';
wSmall = 4; %small scale gaussian filter
wLarge = 8; %large scale gaussian filter
Small = 30; %definition of small objects to remove using bwareaopen
Tau = 0.05; %percentage of tissue area required for analysis (contiguous piece)
MuGlass = 20; %initial guess at glass mean
lMuGlass = 0; %lower limit on glass mean
uMuGlass = 40; %upper limit on glass mean
MuTissue = 110; %initial guess at background tissue mean
lMuTissue = 80; %lower limit on background tissue mean
uMuTissue = 140; %upper limit on background tissue mean
MuNuclei = 200; %initial guess at nuclei tissue mean
lMuNuclei = 160; %lower limit on nuclei mean
uMuNuclei = 240; %upper limit on nuclei mean
SigSqGlass = 4;
SigSqTissue = 10;
SigSqNuclei = 15;
Bins = 0:255; %histogram bins
Opts = optimset('Display','off');

%add paths
addpath '/home/lcoop22/ImageAnalysis/OpenSlide/'
addpath '/home/sgurban/ImageAnalysis/ColorDeconvolution/'
addpath '/home/lcoop22/ImageAnalysis/AdaptiveHMinima/'
addpath '/home/sgurban/Bioinformatics/MixtureModeling/'
addpath '/home/lcoop22/ImageAnalysis/IHCPipeline/'
addpath '/home/lcoop22/ImageAnalysis/PoissonBinarization/'
addpath '/home/lcoop22/ImageAnalysis/MaxFlowBinarization/'
addpath '/home/lcoop22/ImageAnalysis/BoundaryValidator/'
addpath '/home/lcoop22/ImageAnalysis/Utilities/'
addpath '/home/lcoop22/ImageAnalysis/ConstrainedLOG/'

%read in stain calibration file
load(Mfile);

%start timer
tic;

%extract slide name
Dots = strfind(Slide, '.');
Slashes = strfind(Slide, '/');
SlideName = Slide(Slashes(end)+1:Dots(end)-1);

%update console
fprintf('Building statistical models of glass/tissue/nuclei... ');

%check if slide can be opened
Valid = openslide_can_open(Slide);

%slide is a valid file - calculate entropy at 5X to order tiles
if(Valid)

    %get schedule for 5X tiles
    [Level, Scale, Tout, Factor, X, Y, dX, dY] = ...
               TilingSchedule(Slide, ModelingMag, T/(Desired/ModelingMag));

    %read in tiles, calculate entropy to order tiles in terms of content
    tic
    I5 = openslide_read_regions(Slide, repmat(Level, size(X)),...
                        X, Y, repmat(Tout, size(X)), repmat(Tout, size(X)));
    I5 = cellfun(@(x)rgb2gray(x), I5, 'UniformOutput', false);
    Ent5 = cellfun(@(x)entropy(double(x)/255), I5);
    [~, EntOrder] = sort(-Ent5);

    %free up space
    clear I5;

end

% %sample glass/background pixels from designated background tile
% load(BackgroundFile);
% Background = ColorDeconvolution(double(Background), M, [1 1 0], false);
% Scaled = uint8(Background);
% Scaled = imfilter(imcomplement(Scaled(:,:,1)),...
%                 fspecial('gaussian', wLarge, (wLarge/2)/3));
% MuGlass = mean(double(Scaled(BackgroundMask)));
% SigSqGlass = var(double(Scaled(BackgroundMask)));
%
% %generate glass histogram
% HGlass = histc(Scaled(BackgroundMask).', Bins);

%initialize pixel histogram
H = zeros(size(Bins));

%generate schedule for desired magnification
[Level, Scale, Tout, Factor, X, Y, dX, dY] = ...
               TilingSchedule(Slide, Desired, T);

%sample tissue/nuclei pixels from tiles with high entropies
for i = 1:min(5, length(X))

    %read in region corresponding to i'th highest entropy tile
    I = openslide_read_regions(Slide, Level, X(EntOrder(i)), Y(EntOrder(i)), T, T);
    I = I{1};

    %detect large blocks of background value ([0 0 0]) on tile boundary
    Background = (I(:,:,1) == 0) & (I(:,:,2) == 0) & (I(:,:,3) == 0);
    Edge = imclearborder(Background);
    Background = xor(Edge, Background);

    %separate hematoxylin and DAB channels
    Intensity = ColorDeconvolution(double(I), M, [1 1 0], false);
    Scaled = uint8(Intensity);
    Scaled = imfilter(imcomplement(Scaled(:,:,1)),...
                fspecial('gaussian', wLarge, (wLarge/2)/3));

    %make a mask for the tissue region
    Mask = Scaled >= MuGlass + 2*SigSqGlass;
    Mask = imfill(Mask, 'holes');
    Mask = imdilate(Mask, strel('disk', 10));

    %generate histograms for tissue-rich tile
    H = H + histc(Scaled(~Background & Mask).', Bins); %generate histogram

end

%normalize H
H = H / sum(H);

%find peaks in tissue, nuclei limits
TissuePeak = max(H((Bins >= lMuTissue) & (Bins <= uMuTissue)));
NucleiPeak = max(H((Bins >= lMuNuclei) & (Bins <= uMuNuclei)));
GlassPeak = max(H(Bins < lMuTissue));

%optimize model of glass, tissue and nuclei
x0 = zeros(1,8);
x0(1) = GlassPeak / (GlassPeak + TissuePeak + NucleiPeak); %prob(glass) mixture modeling
x0(2) = TissuePeak / (GlassPeak + TissuePeak + NucleiPeak); %prob(tissue) mixture modeling
x0(3) = MuGlass; x0(4) = SigSqGlass;
x0(5) = MuTissue; x0(6) = SigSqTissue;
x0(7) = MuNuclei; x0(8) = SigSqNuclei;
fun = @(x,xdata)x(1)*normpdf(xdata, x(3), x(4)) + ...
        x(2)*normpdf(xdata, x(5), x(6)) + ...
        (1-x(1)-x(2))*normpdf(xdata, x(7), x(8));
GlobalModel = lsqcurvefit(fun, x0, Bins, H,...
                        [0 0 lMuGlass 0 lMuTissue 0 lMuNuclei 0],...
                        [1 1 uMuGlass 10 uMuNuclei 50 uMuNuclei 50], Opts);

%extract parameters from global model
gProbGlass = GlobalModel(1);
gProbTissue = GlobalModel(2);
gProbNuclei = 1-gProbGlass-gProbTissue;
gMuGlass = GlobalModel(3);
gSigSqGlass = GlobalModel(4);
gMuTissue = GlobalModel(5);
gSigSqTissue = GlobalModel(6);
gMuNuclei = GlobalModel(7);
gSigSqNuclei = GlobalModel(8);

% %plot histograms and fit of global model
% figure; plot(Bins, H, 'b'); hold on;
% plot(Bins, fun(GlobalModel,Bins), 'c:');
% plot(Bins, gProbGlass*normpdf(Bins, gMuGlass, gSigSqGlass), 'g');
% plot(Bins, gProbTissue*normpdf(Bins, gMuTissue, gSigSqTissue), 'r');
% plot(Bins, gProbNuclei*normpdf(Bins, gMuNuclei, gSigSqNuclei), 'm');
% xlabel('Distance'); ylabel('frequency');
% xlim([min(Bins) max(Bins)]);
% title('Glass/Tissue/Nuclei');

%define ML transitions for glass / ~glass - test for presence of tissue
MLglass = max(Bins(gProbGlass*normpdf(Bins, gMuGlass, gSigSqGlass) > ...
              gProbTissue*normpdf(Bins, gMuTissue, gSigSqTissue)));
if(isempty(MLglass)) %optimization failed - use lowest point between tissue/glass to threshold
    fprintf('\tOptimization failed - using a-priori model of glass/tissue boundary');
    MLglass = 30;
end

%update console
fprintf('%g seconds.\n', toc);
fprintf('Total tiles: %d\n', length(X));

%process each tile
for i = 1:length(X)

    %start timer
    tic;

    %update console
    fprintf('Processing image %s, X=%s, Y=%s tile %d of %d... ',...
                Slide, num2str(dX(i)), num2str(dY(i)), i, length(X));

    %read in tile
    I = openslide_read_regions(Slide, Level, X(i), Y(i), T, T);
    I = I{1};
%detect large blocks of background value ([0 0 0]) on tile boundary
    Background = (I(:,:,1) == 0) & (I(:,:,2) == 0) & (I(:,:,3) == 0);
    Edge = imclearborder(Background);
    Background = xor(Edge, Background);

    %separate hematoxylin and DAB channels
    Intensity = ColorDeconvolution(double(I), M, [1 1 0], false);
    Scaled = uint8(Intensity);
    Smoothed = imfilter(imcomplement(Scaled(:,:,1)),...         %apply gaussian smoothing
                            fspecial('gaussian', wLarge, (wLarge/2)/3));
    Mask = Smoothed > MLglass;
    Mask(Mask & Background) = false;
    Area = regionprops(Mask, 'Area');
    Max = max([Area.Area]);
    Percentage = Max / sum(~Background(:));
    if(Percentage >= Tau)

        %identify tissue area using global model
        pGlass = gProbGlass*normpdf(Bins, gMuGlass, gSigSqGlass);
        pGlass = pGlass / sum(pGlass);addpath('/home/lcoop22/ImageAnalysis/OpenSlide/');
        pTissue = gProbTissue*normpdf(Bins, gMuTissue, gSigSqTissue) + ...
                    gProbNuclei*normpdf(Bins, gMuNuclei, gSigSqNuclei);
        pTissue = pTissue / sum(pTissue);
        Tissue = MaxFlowBinarization(Smoothed, pGlass, pTissue, 5, 4);
        Tissue = Tissue == 1 & ~Background;

        %tile-centric mixture modeling to separate tissue background and nuclei
        Smoothed = imfilter(imcomplement(Scaled(:,:,1)),...
            fspecial('gaussian', wSmall, (wSmall/2)/3));
        H = histc(Smoothed(Tissue), Bins).'; %generate histogram for tissue area
        H = H / sum(H);
        x0 = zeros(1,5);
        x0(1) = H(round(gMuTissue)) / (H(round(gMuTissue)) + H(round(gMuNuclei)));
        x0(2) = gMuTissue;
        x0(3) = gSigSqTissue;
        x0(4) = gMuNuclei;
        x0(5) = gSigSqNuclei;
        fun = @(x,xdata)x(1)*normpdf(xdata, x(2), x(3)) + (1-x(1))*normpdf(xdata, x(4), x(5));
        TileModel = lsqcurvefit(fun, x0, Bins, H,...
            [0 0.9*gMuTissue 0.9*gSigSqTissue 0.9*gMuNuclei 0.9*gSigSqNuclei],...
            [1 1.1*gMuTissue 1.1*gSigSqTissue 1.1*gMuNuclei 1.1*gSigSqNuclei], Opts);

        %extract parameters from tile-centric model
        tProbTissue = TileModel(1);
        tProbNuclei = 1-tProbTissue;
        tMuTissue = TileModel(2);
        tSigSqTissue = TileModel(3);
        tMuNuclei = TileModel(4);
        tSigSqNuclei = TileModel(5);

%define ML transitions between tissue background and nuclei
        pTissue = tProbTissue*normpdf(Bins, tMuTissue, tSigSqTissue);
        pNuclei = tProbNuclei*normpdf(Bins, tMuNuclei, tSigSqNuclei);
        Nuclei = MaxFlowBinarization(Smoothed, pTissue, pNuclei, 4, 4);
        Nuclei = Nuclei == 1;
        Nuclei = Nuclei & Tissue;

        %perform multiscale log filter to enhance nuclear regions
        F1 = MultiscaleLog(Scaled(:,:,1), 5, 30); %apply multiscale LOG filter
        F1 = F1 / max(F1(:)); %scale to interval (0, 1)
        Nuclei = Nuclei & (F1 >= 0.1);
        Nuclei = bwareaopen(Nuclei, Small, 4);
        [Marker, hadp] = AdaptiveHMinima2(Nuclei, 4);
        Distance = bwdistgeodesic(Nuclei, Marker);
        Distance(isnan(Distance)) = 0;
        L = watershed(Distance);

        %eliminate small objects and quantify stain in each object
        Props = regionprops(L, Scaled(:,:,2), 'Area');
        [~, BG] = max([Props.Area]);
        L(L == BG) = 0; Props(BG) = [];
        L(L > BG) = L(L > BG) - 1;
        Filter = [Props.Area] < Small;
        L(ismember(L, find(Filter))) = 0;

        %check if objects were discovered
        if(max(L(:)) > 0)

            %capture boundary info
            [Bounds, L] = bwboundaries(L>0, 8, 'noholes');
            bX = cellfun(@(x)x(:,2), Bounds, 'UniformOutput', false);
            bY = cellfun(@(x)x(:,1), Bounds, 'UniformOutput', false);

            %capture object centroids, features
            Props = regionprops(L, Scaled(:,:,2), 'Area', 'MeanIntensity',...
                'MinIntensity', 'MaxIntensity', 'Centroid');
            Centroids = cat(1,Props.Centroid);
            cX = Centroids(:,1);
            cY = Centroids(:,2);
            Features = [[Props.MeanIntensity]; [Props.MinIntensity]; [Props.MaxIntensity]];

            %add scan and analysis magnifications, generate names
            ScanMag = repmat(Desired/Scale, [1 size(Features,2)]);
            AnalysisMag = repmat(Desired, [1 size(Features, 2)]);
            Features = [ScanMag; AnalysisMag; Features];

            %Saumya - add area
            Features = [Features; [Props.Area]];
            FeatureNames = {'ScanMag', 'AnalysisMag', 'Mean.DAB.Intensity', 'Min.DAB.Intensity', 'Max.DAB.Intensity', 'Cell Area (px)'};
            
            %place boundaries, centroids in global frame
            cX = cX + dX(i);
            cY= cY + dY(i);
            bX = cellfun(@(x) x + dX(i), bX, 'UniformOutput', false);
            bY = cellfun(@(x) x + dY(i), bY, 'UniformOutput', false);

            %merge colinear points on boundaries
            for k = 1:length(bX)
                [bX{k}, bY{k}] = MergeColinear(bX{k}, bY{k});
            end

            %parse filename for generating output (remove path, file extension)
            Slashes = strfind(Slide, '/');
            Dots = strfind(Slide, '.');
            if(isempty(Slashes))
                Slashes = 0;
            end
            Xstr = sprintf('%06.0f', dX(i));
            Ystr = sprintf('%06.0f', dY(i));

            %save features, generate database text file, and visualization
            save([Folder SlideName '.' Xstr '.' Ystr '.mat'],...
                    'L', 'Features', 'cX', 'cY');
            SegmentationReport([Folder SlideName '.seg.' Xstr '.' Ystr '.txt'],...
                    SlideName, cX, cY, double(Features.'), FeatureNames, bX, bY);
            imwrite(EmbedBoundaries(I, LabelPerim(L, 4), [0 1 1]),...
                [Folder SlideName '.seg.' Xstr '.' Ystr '.jpg']);

        end

    end

    %update console
    fprintf('%g seconds.\n', toc);

end
