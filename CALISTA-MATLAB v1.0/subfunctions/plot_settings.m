function plot_settings(xvec,yvec,xLABEL,yLABEL)
ax = gca;
set(ax,'ytick',yvec)
set(ax,'yticklabel',yLABEL)%(1:tickStep:numel(genes)))

set(ax,'xtick',xvec)
set(ax,'xticklabel',xLABEL)%(1:tickStep:numel(genes)))
ax.XTickLabelRotation = 45; %xlabel rotation
set(ax, 'Ticklength', [0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dx=xvec(2)-xvec(1);
   dy=yvec(2)-yvec(1);

   nx=numel(xvec)+1; 
   ny=numel(yvec)+1; 

   
   x = linspace(min(xvec)-dx/2,max(xvec)+dx/2,nx);
   x = [repmat(x',1,ny)'; nan(1,nx)];
   xx=x(1:ny,:)';
   xx=[xx ; nan(1,ny)];
   y=[linspace(min(yvec)-dy/2,max(yvec)+dy/2,ny) nan]';
   y = repmat(y,1,nx);
   yy=y(1:ny,:)';
   yy=[yy; nan(1,ny)];

   
   xx=xx(:);
   yy=yy(:);
   x = x(:);
   y = y(:);
   
   hold on;
   line(x,y,'linestyle','-','color','k');
   line(xx,yy,'linestyle','-','color','k')
%%%%%%%%%%%%%%%%%%%%%

set(ax, 'layer', 'top');
end