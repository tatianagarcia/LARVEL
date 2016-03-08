% Preallocate the array
I=imread(fullfile(pathname,filename{1}));
Images=zeros([size(I) nFrames],class(I));%prealocate memory
Images(:,:,1)=rgb2gray(I);
% Create image sequence array
for p=2:nFrames
   Images(:,:,p)=rgb2gray(imread(fullfile(pathname,filename{p})));
end
clear p filename pathname