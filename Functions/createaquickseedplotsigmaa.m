function [h] = createaquickseedplotsigmaa(iIndex,truesigmaa,outsl,outfull_sw,outfull_psd_estrank,outfull_psd_truerank,mypos,face,x,y,z,depthLabels)

h=subplot(3,3,1);
trisurf(face,x,y,z,truesigmaa(:,iIndex),'EdgeColor','none')
%colorbar()
title('Truth')
colorrange=caxis();
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;

subplot(3,3,2)
h=trisurf(face,x,y,z,truesigmaa(:,iIndex),'EdgeColor','none');
caxis(colorrange);
colorbar();
h.FaceAlpha=0;
axis off;
axis equal;

h=subplot(3,3,4);
trisurf(face,x,y,z,outsl.smSA_symm(:,iIndex),'EdgeColor','none')
%colorbar()
caxis(colorrange);
title(depthLabels{3})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,5);
trisurf(face,x,y,z,outsl.smSA_psd(:,iIndex),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{4})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,6);
trisurf(face,x,y,z,outfull_sw.smSA_symm(:,iIndex),'EdgeColor','none')
caxis(colorrange)
title(depthLabels{1})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,7);
trisurf(face,x,y,z,outfull_sw.smSA_psd(:,iIndex),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{2})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,8);
trisurf(face,x,y,z,outfull_psd_truerank.smSA_psd(:,iIndex),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{6})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;


h=subplot(3,3,9);
trisurf(face,x,y,z,outfull_psd_estrank.smSA_psd(:,iIndex),'EdgeColor','none')
caxis(colorrange);
title(depthLabels{5})
axis equal;
axis off;
a = text(x(iIndex),y(iIndex),z(iIndex), 'o','HorizontalAlignment','center','VerticalAlignment','middle');
a.FontSize=15;
h.CameraPosition = mypos;

end

