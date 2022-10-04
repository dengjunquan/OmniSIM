
% -----------------------------------------------------
% -- Fast mmWave Ray Tracing Simulator (v0.2)
% -- 2018 (c) junquan.deng@aalto.fi
% -----------------------------------------------------
close;clear;clc;
rng(888);
%%
filename_walls = 'data/Buildings.json';
filename_roads = 'data/UEs.json';
filename_BSs   = 'data/BS.json';
filename_Scatters = 'data/scatters.json';

[walls,polygons,centerLonLat] = FindWalls(filename_walls);
roads                         = FindRoads(filename_roads,centerLonLat);
BSposAll                      = FindBSs(filename_BSs,centerLonLat);
Scatters                      = FindBSs(filename_Scatters,centerLonLat);

UEdensity = 2;
[UEs, TimeStamps]             = FindUEsonRoads(roads,UEdensity);
UEs = UEs';
BSOrientationAll = ...        % BS Orientations
    [0
    ];

wallsMatrix = zeros(length(walls),4);
wallsCenter = zeros(length(walls),2);

for walli = 1:length(walls)
    wallsMatrix(walli,1) = walls{walli,1}.p1(1);
    wallsMatrix(walli,2) = walls{walli,1}.p1(2);
    wallsMatrix(walli,3) = walls{walli,1}.p2(1);
    wallsMatrix(walli,4) = walls{walli,1}.p2(2);
    wallsCenter(walli,1) = (wallsMatrix(walli,1) + wallsMatrix(walli,3))/2;
    wallsCenter(walli,2) = (wallsMatrix(walli,2) + wallsMatrix(walli,3))/2;
end

ToPlot = 1;

if ToPlot ==1
    figure; hold on;
    % Draw walls
    if ~isempty(walls)
        for i=1:length(walls)
            plot([walls{i}.p1(1),walls{i}.p2(1)],[walls{i}.p1(2),walls{i}.p2(2)],'b');
        end
    end
    axis tight; axis equal;
    grid on;
end

xlim([-50  350]);
ylim([-250 150]);

if ToPlot == 1
    % Draw the UEs
    scatter(UEs(:,1), UEs(:,2), '.');
end

% Draw BSs
if ToPlot ==1
    for bi=1:size(BSposAll,1)
        BSxy(1) = BSposAll(bi,1);
        BSxy(2) = BSposAll(bi,2);
        plot(BSxy(1),BSxy(2),'r^','linewidth',1); hold on;
        %text(BSxy(1) + 10,BSxy(2),num2str(bi),'FontSize',12);
    end
end

% Draw Scatters
if ToPlot ==1
    for si=1:size(Scatters,1)
        plot(Scatters(si,1),Scatters(si,2),'g*','linewidth',1); hold on;
        %text(BSxy(1) + 10,BSxy(2),num2str(bi),'FontSize',12);
    end
end

bsindex = 1;
BSxy(1) =  BSposAll(bsindex,1);
BSxy(2) =  BSposAll(bsindex,2);

%%
BouncingOrder = 4;
UEindex = 50;
testUEs = 1:length(UEs);
numBSs = size(BSposAll,1);
path_arrays_allBSs  = cell(numBSs);

for BSindex = 1:numBSs
    disp(['BS = ',num2str(BSindex)]);
    BSpos = BSposAll(BSindex,:);
    path_arrays = FastRT(walls,BSpos,Scatters,UEs,BouncingOrder);  % Ray tracing based on shooting-and-bouncing
    path_arrays_allBSs{BSindex} = path_arrays;
end
%%
numUEs = length(UEs);

index = 10;
path_array = path_arrays{index};

for pathi = 1:length(path_array)
    path  = path_array{pathi};
    pathlen = size(path,1);
    for segmenti = 1:pathlen-1
        xy1 = path(segmenti,1:2);
        xy2 = path(segmenti+1,1:2);
        plot([xy1(1) xy2(1)],[xy1(2) xy2(2)],'-','Color',[0.5 0.5 0.5]);
    end
end

plot(UEs(UEindex,1), UEs(UEindex,2), 'rO');
plot(BSxy(1),BSxy(2),'b^'); hold on;

%[~,SortedIndexs] = sort(UEs(:,2));  % for coloring purpose
%UEs = UEs(SortedIndexs,:);
%TimeStamps = TimeStamps(SortedIndexs);
%% ======================== Channel realization ===============================
BSobj.ElementPattern = 'patchcosine'; % 1- 'omnidirectional'; 2 - 'patchcosine'; 3 - 'sector';
BSobj.array          = 'UPA';

M=16; N=1;  % equivalent to a ULA

if strcmp(BSobj.array,'UPA')==1
    BSobj.M = M;
    BSobj.N = N;
    BSobj.dx = 1/2;
    BSobj.dz = 1/2;
end

UEobj.array = 'UPA';
M=1;
N=1;

UEobj.ElementPattern = 'omnidirectional';
if strcmp(UEobj.array,'UPA')==1
    UEobj.M = M;
    UEobj.N = N;
    UEobj.dx = 1/2;
    UEobj.dz = 1/2;
end

numRealisation = BSobj.M*2;   % number of channel realizations

Pbs = 10^(33/10);      % TX power
Bandwidth   = 1e8;
NoiseFloor  = -174;
NoiseFigure = 6;
NoisePower  = Bandwidth*10^((NoiseFloor + NoiseFigure)/10);

numWalls = length(walls);
eta_rs   = unifrnd(4,6,1,numWalls); % random relative permittivity of walls;

Nsubrays    = 20;         % # of sub-rays in each cluster, we use one sub-ray here
sigma_phy   = 5/180*pi;  % standard deviation of the angles in elevation both of Rx and Tx
sigma_theta = 10/180*pi;  % standard deviation of the angles in azimuth   both of Rx and Tx
max_theta   = 2*sigma_theta;
max_phy     = 2*sigma_phy;

%==========================================================================
numtestUEs   = length(testUEs);
% The channel realisation is stored in Hs_t
Hs_t         = zeros(numRealisation, numBSs, numtestUEs, BSobj.M*BSobj.N, UEobj.M*UEobj.N);
UEsforEachBS = cell(numBSs,1);
Noise = sqrt(NoisePower/Pbs)*sqrt(1/2)*(randn(1,BSobj.M)+1j*randn(1,BSobj.M));  % Gaussian noise
NoiseLevel = norm(Noise);

ChObjectsAll = cell(numBSs,1);
for BSindex = 1:numBSs
    Orientation = BSOrientationAll(BSindex);
    path_arrays = path_arrays_allBSs{BSindex};
    for iter = 1:numRealisation
        disp(['BS =',num2str(BSindex),', Channel Realisation = ', num2str(iter)]);
        if iter ==1
            [Hs,ChObjects]   = SalehValenzuela(path_arrays(testUEs),walls,eta_rs,sigma_theta,sigma_phy,max_theta,max_phy,Nsubrays,BSobj,UEobj,Orientation);
            UEsforEachBS{BSindex}  = testUEs(mean(abs(Hs),2) > NoiseLevel/2);  % Only UEs with sufficient receiverd powers are heard by the BS
        else
            Hs  = ChannelUpdate(ChObjects,numtestUEs,Nsubrays,BSobj,UEobj);    % Update of fast fading
        end
        Hs_t(iter, BSindex, :,:,:)  = Hs;
    end
    ChObjectsAll{BSindex} = ChObjects;
end

%%
ChObjects = ChObjectsAll{1};
Object  = ChObjects{10};
N = length(Object);
figure; hold on;
for raynow = 1:N
    %if Object{raynow}.LOS==0
        subpathdelays = Object{raynow}.delay + Object{raynow}.subrayDelays;
        subpathpowers = Object{raynow}.ClusterPower.*Object{raynow}.subrayPowers;
        subpathDoAs   = Object{raynow}.DoA + Object{raynow}.DoAs_theta;
     
        x = subpathdelays;
        y = subpathDoAs;
        z = 10.*log10(subpathpowers);
        scatter3(x,y,z,'ro','filled');
        zlim([-180,-50]);
        view(-30,30); grid on;
    %end
end
for raynow = 1:N
    %if Object{raynow}.LOS==0
        subpathdelays = Object{raynow}.delay + Object{raynow}.subrayDelays;
        subpathpowers = Object{raynow}.ClusterPower.*Object{raynow}.subrayPowers;
        subpathDoAs   = Object{raynow}.DoA + Object{raynow}.DoAs_theta;
        x = subpathdelays;
        y = subpathDoAs;
        z = 10.*log10(subpathpowers);
        zlim([-180,-50]);
        zRng = zlim;
        for k = 1:length(x)
            xL = [x(k) x(k)];
            yL = [y(k) y(k)];
            zL = [zRng(1) z(k)];
            plot3(xL,yL,zL,'k-','linewidth',0.25)         % plot vertical line (3D)
        end
        view(-30,30); grid on;
    %end
end

xlabel('Delay');
ylabel('Azimuth DoA');
zlabel('Path Power');
%%
channelpowers = zeros(1,numtestUEs);
for ui=1:numtestUEs
    hnorm = norm(squeeze(Hs_t(1,1,ui,:,:)));
    channelpowers(ui) = 10.*log10(hnorm^2);
end
figure;
plot(1:numtestUEs,channelpowers,'linewidth',2);
grid on;
xlabel('Time in second');
ylabel('Channel power $|\mathbf{H}|^2$ in dB','interpreter','latex');