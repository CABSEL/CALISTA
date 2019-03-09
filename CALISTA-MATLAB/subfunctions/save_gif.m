function save_gif(fig,filename,num_subplots,p1,p2)
% save gif
figure(fig)
az = 0;
el = 45;
view([az,el])
degStep = 5;
detlaT = 0.00001;
fCount = 71;
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,fCount) = 0;
k = 1;
del=0.1;
% spin 45°
count=1;
%filename='aaa.gif';
for ii = 0:2:360
    az = ii;
    for i=1:num_subplots
        subplot(p1,p2,i)
        view([az,el])
    end
        
    frame = getframe(fig); %or gcf?
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if count == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
    end
    count=count+1;
end
%imwrite(im,map,filename2,'DelayTime',0,'LoopCount',inf)
