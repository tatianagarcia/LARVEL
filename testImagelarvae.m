close all
clear all
% P9160004.jpg
% using examples from:
% https://www.mathworks.com/help/images/noise-removal.html
% https://www.mathworks.com/help/imaq/examples/calculating-the-length-of-a-pendulum-in-motion.html

% load the image
a=imread('DSC_0337.jpg');

figure();
imagesc(a);  %,[2^0 2^12])
colormap gray;
%%
% pass from RGB to gray image 
b1=rgb2gray(a);
figure();
imagesc(b1);  %,[2^0 2^12])
colormap gray;
%%
% trim the edges (define bottom and free surface, once per series)
b=b1(:,630:4000);
figure();
imagesc(b);   %,[2^0 2^12])
colormap gray;
%%
% original = b in intensity

% now shift colors
c=imcomplement(b);
figure();
imagesc(c,[100 182]);
colormap gray;
% %%
% min(min(c1))
% max(max(c1))
% %%
% c2=2*(c1-min(min(c1)));
% figure()
% imagesc(c)
% colormap gray
% %%
% c3=2*(c2-min(min(c2)));
% c4=2*(c3-min(min(c3)));
% figure()
% imagesc(c2)
% colormap gray
% min(min(c))
% max(max(c))
%%
% now remove specks:

d=wiener2(c,[10 10]); %2D adaptive noise-removal filtering
figure()
imagesc(d,[100,200])
colormap gray
% remove particles (i.e., obtain background)
% note: we should have a picture without particles to subtract first, for
% now we are subtracting a median filtered image
%%
e=medfilt2(d,[50 50]); % 40 is larger than particle diamater
figure()
imagesc(e)
colormap gray
%%
figure()
subplot(121)
imagesc(d,[1 200])
colormap gray

subplot(122)
imagesc(e,[1 200])
colormap gray

%%  get image with only particles

f1=1*(d-e);
figure()
imagesc(f1,[1 10])
colormap gray

% more cleaning

f=wiener2(f1,[40 40]);
figure()
imagesc(f,[1 10])%,[100,200])
colormap gray
%%
% f=d;
% change to binary (1 or 0)

% close all

level = graythresh(f);
g = im2bw(f,0.01); %  0.01 will change based on the imaqe quality
figure()
imagesc(g)%,[1 10])
colormap gray
%%


% more cleaning of non-particles
h=imclearborder(g);
structDisk = strel('disk', 6); % 4 to be adjusted
i = imopen(h, structDisk);

% and find centroids

gcent=regionprops(i,'centroid');
gcentroids=cat(1,gcent.Centroid);

figure()
subplot(241)
imshow(b)
title('original')
subplot(242)
imshow(c)
title('inverted')
subplot(243)
imshow(d)
title('without noise')
subplot(244)
imshow(e)
title('background -particles removed')
subplot(245)
imshow(f)
title('only particles')
subplot(246)
imshow(g)
title('bw')
subplot(247)
imshow(i)
title('bw cleaned')
subplot(248)
imshow(b)
hold on
plot(gcentroids(:,1),gcentroids(:,2), 'r*')
hold off
title('centroids')

figure()
imshow(b)
hold on
plot(gcentroids(:,1),gcentroids(:,2), 'r*')
hold off
title('centroids')