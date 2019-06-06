function [] = EventsSurfacePlot( Lall,event )

writerObj = VideoWriter('events.avi');
    writerObj.FrameRate = 1;
    open(writerObj);
    for k = 1:min(10,event)
	SurfacePlot(Lall,k)
    frame=getframe(gcf);
    writeVideo(writerObj,frame);
    end
    close(writerObj);
    
end