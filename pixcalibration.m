function [Pix_Cmx,Pix_Cmy,Xpixels,Ypixels,Xcm,Ycm]=pixcalibration(I)
n=input('number of points to use in the calibration ');
%n_sq_x=3; n_sq_y=4;
%n=20
dx=1;%1cm
dy=1;%1cm
h=figure(1);
imshow(I,'InitialMagnification','fit');
set(h,'color',[1 1 1],'position',[ 1 1 1920 1080]);
Annotation=annotation(h,'textbox',[0.45 0.8 0.28 0],...
    'String',{'Zoom in and click in consecutive points'},...
    'FontWeight','bold','FontSize',12,'FitBoxToText','off',...
    'EdgeColor','none','Color',[1 1 1]);
title('Click on the figure every 1cm');
xlabel('X(Pixels)');ylabel('Y(Pixels)');
disp('Click on the figure every 1cm');
x= [];
y = [];
figure(1);
hold on;
pause
for count=1:n+1,
    if count==n+1
        delete(Annotation)
        Annotation=annotation(h,'textbox',[0.45 0.8 0.28 0],...
            'String',{'Click in a reference point of the grid'},...
            'FontWeight','bold','FontSize',12,'FitBoxToText','off',...
            'EdgeColor','none','Color',[1 1 1]);
        title('Locate a reference point in grid','color','r');
        [Xpixels,Ypixels] = ginput(1);
        plot(Xpixels,Ypixels,'o','color','r','linewidth',3);
        pause(1)
        title('Please go to the command window','color','k','fontsize',16);
        %close(h)
        %ginput4???
    else
        [xi,yi] = ginput(1);
        figure(1);
        plot(xi,yi,'+','color',[ 1.000 0.314 0.510 ],'linewidth',2);
        x = [x;xi];
        y = [y;yi];
    end
    
end;
%pixPerCm=Diff(y)/1;
Xcm=[];Ycm=[];
while isempty(Xcm)
clc
Xcm=input('Coordinate X in cm: ');
end
while isempty(Ycm)
clc
Ycm=input('Coordinate Y in cm: ');
end
pixels=[x-min(x) y-min(y)];
clear xi yi count x y
%% Consequtive Points in a grid e.g.
% * * *
% * * *
n_points_x_default = 4; n_points_y_default = 5;
n_sq_x = input(['Number of points (min 2) along the X direction ([]=' num2str(n_points_x_default) ') = '])-1; %6
if isempty(n_sq_x), n_sq_x = n_points_x_default-1; end;
n_sq_y = input(['Number of points (min 1) along the Y direction ([]=' num2str(n_points_y_default) ') = '])-1; %6
if isempty(n_sq_y), n_sq_y = n_points_y_default-1; end;
% Points in a ruler with constant x=0
if n_sq_y==0
    cm=[(0:dy:n-1)' zeros(n,1)];
else
    y_l = ((0:n_sq_y)'*ones(1,n_sq_x+1)); x_l = (ones(n_sq_y+1,1)*(0:n_sq_x));
    cm= [x_l(:) y_l(:)];
    if (n_sq_y+1)*(n_sq_x+1)>n
        cm(n+1:(n_sq_y+1)*(n_sq_x+1),:)=[];
    end
end
%Least squares pixels=cm*A --> A=[pixeles/cm]
A=lscov(cm,pixels);%To get cm-->cmx=pixels(:,1)/A(1,1);
Pix_Cmx=A(1,1);
Pix_Cmy=A(1,2);
close(h)


