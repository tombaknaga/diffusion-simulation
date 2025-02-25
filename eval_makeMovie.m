% load init data
load("initConfigData.mat");

doVideo = false;

if doVideo
    % initialization for video output
    videofile = 'diffusionMovie.mp4';
    videoWriter = VideoWriter(videofile, 'MPEG-4');
    videoWriter.FrameRate = 60;
    open(videoWriter);
end

drawtime = 0;

files = dir("dat/*.mat");
files = natsortfiles({files.name});

resX = 1280;
resY = 720;


for i = 1:length(files)
    filename = horzcat('dat/', files{i});
    load(filename);
    disp(filename);
    % Capture and write frame
    if t(i) >= drawtime

        figure(123)
        clf
        set(gca,'looseinset',get(gca,'tightinset'),'fontsize',16)
        set(gcf,'position',[0 0 resX resY])
        hold on
        surf(X, Y, u, EdgeColor='none');
        contour3(X, Y, u,25,Edgecolor='k',EdgeAlpha=0.40);
        hold off
        title(['Time: ', t, ' s']);
        xlabel('X [m]'); ylabel('Y [m]'); zlabel('Concentration');
        zlim([0 maxConMult])
        pbaspect([1 1 0.2])

        clim([0 maxConMult]);

        view(3);

        drawnow;
        drawtime = drawtime + 1/60;

        if doVideo
            frame = getframe(gcf);
            writeVideo(videoWriter, frame);
            
        end
        
    end
end

% Close Video Writer

if doVideo
    close(videoWriter);
    disp(['Simulation video saved as ', videofile]);
end

