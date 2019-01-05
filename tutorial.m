%% TUTORIAL: Lena Analysis 



%% load Lena image 

load('lena512.mat')
I = uint8(lena512);
figure;imagesc(I);colormap('gray');axis off



%% extract the (horizontal) image visibility graph (IHVG)

Edge_list=imageVisibilityGraph(I,'horizontal',true);

G = graph(Edge_list(:,1),Edge_list(:,2));

figure; plot(G);axis off;title('Lena IHVG') 

% we could neglect the regular lattice structure that all IVGs have by construction
% by setting the <lattice> parameter equal to false:
% Edge_list=imageVisibilityGraph(I,'horizontal',false);

%% apply IHVG k-filter (degree filter) to the Lena image

% degree sequence from the IHVG
Deg_seq = degree(G);

% Degree distribution P(k)
Pk=hist(Deg_seq,1:1:max(Deg_seq))./(size(I,1)*size(I,2));

% Pixel intensity distribution from original image
Pi=hist(double(I(:)),1:1:max(double(I(:))))./(size(I,1)*size(I,2));


%filter construction
I2 = reshape(Deg_seq,size(I,1),size(I,2));
I2 = uint8(full(I2'));


figure;
subplot(2,2,1);imagesc(I);colormap('gray');axis off;title('Lena original');
subplot(2,2,2);plot(1:1:max(double(I(:))),Pi);xlabel('Pixel Intensity');ylabel('Frequency');
subplot(2,2,3);imagesc(I2);colormap('gray');axis off;title('k-filter');
subplot(2,2,4);plot(1:1:max(Deg_seq),Pk);xlabel('Pixel Intensity');ylabel('Frequency');


%% extract image visibility patches for texture classification


load('texture_examples.mat')


figure;
subplot(2,2,1);imagesc(canvas);colormap('gray');axis off;title('Canvas');
subplot(2,2,2);imagesc(cushion);colormap('gray');axis off;title('Cushion');
subplot(2,2,3);imagesc(linseeds);colormap('gray');axis off;title('Lin Seeds');
subplot(2,2,4);imagesc(sand);colormap('gray');axis off;title('Sand');


%%
% number of pixels/nodes in the graphs 
N=size(canvas,1).^2;

Iseq={canvas,cushion,linseeds,sand};

Pk=cell(1,numel(Iseq));
Z=cell(1,numel(Iseq));

% Set a maximum degree 
maxK=100;

for i=1:numel(Iseq)

Edge_list=imageVisibilityGraph(Iseq{i},'horizontal',true);
Deg_seq = hist(Edge_list(:,1),1:N)+hist(Edge_list(:,2),1:N);

% degree distribution is extracted as a global-image-properties descriptor
Pk{i}=hist(Deg_seq,1:1:maxK)./N;

% patch frequency is extracted as a local-image-properties descriptor 
% we have set stride parameter equal 1
Z{i}=visibilityPatches(Iseq{i},1,'horizontal');

end



figure;

subplot(2,1,1);
loglog(1:1:maxK,Pk{1},'-*');
hold on
loglog(1:1:maxK,Pk{2},'-*');
loglog(1:1:maxK,Pk{3},'-*');
loglog(1:1:maxK,Pk{4},'-*');
xlim([8,maxK])
xlabel('k');ylabel('Frequency');
legend('cushion','canvas','lin seeds','sand')
legend boxoff

subplot(2,1,2);
loglog(Z{1},'-x');
hold on
loglog(Z{2},'-x');
loglog(Z{3},'-x');
loglog(Z{4},'-x');
xlabel('patch id');ylabel('Frequency');
legend('cushion','canvas','lin seeds','sand')
legend boxoff

