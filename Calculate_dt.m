dt=zeros(nFrames,1);%Allocate memory.  numel-->Number of elements in array
for i=1:nFrames
info(i)=imfinfo(fullfile(pathname,filename{i}));%imfinfo->gets the file metadata.  fullfile-->generates the full path for every figure
end
DigCamera=arrayfun(@(x) x.DigitalCamera,info); date={DigCamera.DateTimeOriginal}';
Milliseconds=str2double({DigCamera.SubsecTimeOriginal}')/100;
%str2double-->Converts string to a doble
datevect=datevec(date,'yyyy:mm:dd HH:MM:SS');
datevect(:,end)=datevect(:,end)+Milliseconds;
seconds=(datevect(:,4)*3600)+(datevect(:,5)*60)+datevect(:,end);
dt=diff(seconds);
Time(1)=0;Time(2:nFrames)=cumsum(dt);
