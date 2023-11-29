classdef GridBuildTools

    properties

        PIP = 10 % pitch in pixels = lensPitch/pixelPitch = 1.4000e-05/1.4000e-06
        FDRM = 1/1.5 % Filter Disk Multiply Parameter
        Crop = [25,25] % crop white image in pixels
        Crop_Multi = 2;
        StepSize = 250 % step size during searching peak points
        IM_size; % image size
        FirstPosShiftRow; % Row shift for Microlens
        Rotation; % Lens array rotation angle
        XSpacing; % lens spacing in X
        YSpacing; % lens spacing in Y
        Offset; % offset when building a lens grid model (corrected by adding EstError)
        EstError % estimation error between real grid and constructed grid 
% (only record it and this error has been already added in function PT = ProEst(PT,WhiteImg)
% 'total offset estimation' part
        
    end

    methods 

        % lens grid property estimation including tilting angle, spacing
        % and EstError

        function PT = ProEst(PT,WhiteImg)
            
            filter = fspecial( 'disk', PT.PIP/2  * PT.FDRM );
            filter = filter ./ max( filter(:) ); 
%             imagesc(filter)

            % disk is an odd matrix: conv2 size is the sum of the size of image and the
            % radius of filter (left, right side / up, down side).
            
            Conv2RS_size = size(WhiteImg)+2*(size(filter)-1)/2;

            % FFT of WhiteImg & filter

            WhiteImg_FT = fft2(WhiteImg,Conv2RS_size(1),Conv2RS_size(2));
            filter_FT = fft2(filter,Conv2RS_size(1),Conv2RS_size(2));

            % Conv2RS = conv2 response
            Conv2RS_FT = WhiteImg_FT.*filter_FT;
            Conv2RS = ifft2(Conv2RS_FT);
            CropRegion = (size(Conv2RS)-size(WhiteImg))/2; 
            Conv2RS = Conv2RS(1+CropRegion(1):end-CropRegion(1),1+CropRegion(2):end-CropRegion(2));
            Conv2RS = (Conv2RS  - min(Conv2RS(:))) ./ abs(max(Conv2RS(:))-min(Conv2RS(:)));  
            
            % peak detection & cropping
            Peaks = imregionalmax(Conv2RS);
            PeakId = find(Peaks==1);
            [PeakIdX,PeakIdY] = ind2sub(size(Peaks),PeakId);

            % plot 
%             imagesc(Conv2RS)
%             hold on
%             plot(PeakIdY,PeakIdX,'r.','MarkerSize', 13);
%             set(gca, 'FontSize', 15);
%             axis equal

            PT.IM_size = size(Conv2RS);
            Conv2RS_size = PT.IM_size;
            InsidePts = find(PeakIdX > PT.Crop(1) & ...
                PeakIdX < Conv2RS_size(1) - PT.Crop(1) & ...
                PeakIdY > PT.Crop(2) & ...
                PeakIdY < Conv2RS_size(2) - PT.Crop(2));

            PeakIdX = PeakIdX(InsidePts);
            PeakIdY = PeakIdY(InsidePts);
            
%             imagesc(Conv2RS)
%             hold on
%             plot(PeakIdY,PeakIdX,'r.','MarkerSize', 10);
         
            Triangulation = delaunayTriangulation(PeakIdY, PeakIdX);

            % ********************* Y searching *****************************
            XId = 1;
            
            for X = PT.Crop_Multi*PT.Crop(1) : PT.StepSize : Conv2RS_size(1) - PT.Crop_Multi*PT.Crop(1)
                
                Pos = [X,PT.Crop_Multi*PT.Crop(2)]; % start point for each row
                YId = 1;
            
                while (Pos(2) <= Conv2RS_size(2)-PT.Crop_Multi*PT.Crop(2))
            
                    ClosestLabel = nearestNeighbor(Triangulation, Pos(2), Pos(1));
                    ClosestPt = [PeakIdX(ClosestLabel),PeakIdY(ClosestLabel)];
            
                    RecPtsY(XId,YId,:) = ClosestPt; 
            
                    Pos = ClosestPt+[0,PT.PIP];
            
                    YId = YId+1;
            
                end
                
                LineFitY(XId,:) = polyfit(RecPtsY(XId,1:end-1,2), RecPtsY(XId,1:end-1,1), 1);
            
                XId = XId+1;
                
            end

            % ********************* X searching *****************************            
            YId = 1;
            
            for Y = PT.Crop_Multi*PT.Crop(2) : PT.StepSize : Conv2RS_size(2) - PT.Crop_Multi*PT.Crop(2)
                
                Pos = [PT.Crop_Multi*PT.Crop(1),Y]; % start point for each row
                XId = 1;
            
                while (Pos(1) <= Conv2RS_size(1)-PT.Crop_Multi*PT.Crop(1))
            
                    ClosestLabel = nearestNeighbor(Triangulation, Pos(2), Pos(1));
                    ClosestPt = [PeakIdX(ClosestLabel),PeakIdY(ClosestLabel)];
            
                    RecPtsX(XId,YId,:) = ClosestPt; 
            
                    Pos = ClosestPt+[PT.PIP*sqrt(3),0];
            
                    XId = XId+1;
            
                end
                
                LineFitX(YId,:) = polyfit(RecPtsX(1:end-1,YId,1), RecPtsX(1:end-1,YId,2), 1);
            
                YId = YId+1;
                
            end

            % ********************* estimation *****************************  
            PT.Offset = squeeze(RecPtsY(1,1,:))';
            RecPtsY = RecPtsY(:, 1:end-1,:);  
            RecPtsX = RecPtsX(1:end-1, :,:);
            
            % angle estimation
            AngleX = atan(mean(LineFitX(:,1)));
            AngleY = atan(-mean(LineFitY(:,1)));
            PT.Rotation = mean([AngleX,AngleY]);
            
            % spacing estimation
            X_dist = diff(squeeze(RecPtsX(:,:,1)),1,1);
            X_dist = mean(X_dist(:));

            Y_dist = diff(squeeze(RecPtsY(:,:,2)),1,2);
            Y_dist = mean(Y_dist(:));
          
            PT.XSpacing = 1/2*X_dist / cos(PT.Rotation);
            PT.YSpacing = Y_dist / cos(PT.Rotation);

            % estimation error          
            YMax = floor((Conv2RS_size(2)-PT.Crop_Multi*PT.Crop(2)-PT.Offset(2))/(cos(PT.Rotation)*PT.YSpacing))+1 ;
            XMax = floor((Conv2RS_size(1)-PT.Crop_Multi*PT.Crop(1)-PT.Offset(1))/(cos(PT.Rotation)*PT.XSpacing))+1 ;

            Rot_matrix = [cos(PT.Rotation),-sin(PT.Rotation);sin(PT.Rotation),cos(PT.Rotation)];           
            RotTrans_matrix = [Rot_matrix,PT.Offset'];

            [gridX,gridY] = ndgrid((0:XMax-1).*PT.XSpacing, (0:YMax-1).*PT.YSpacing);
            gridY(2:2:end,:) = gridY(2:2:end,:) + 1/2.*PT.YSpacing;
            GridCoords = (RotTrans_matrix*[gridX(:),gridY(:),ones(numel(gridX),1)]')';
            
            RealPtsId = nearestNeighbor(Triangulation, round(GridCoords(:,2)), round(GridCoords(:,1)));           
            RealPtsCoords = [PeakIdX(RealPtsId), PeakIdY(RealPtsId)];

            PT.EstError  = median(RealPtsCoords - GridCoords);

            % total offset estimation
            Offset_T = [PT.Offset(1)+PT.EstError(1),PT.Offset(2)+PT.EstError(2)];
            InvRot_matrix = [cos(-PT.Rotation),-sin(-PT.Rotation);sin(-PT.Rotation),cos(-PT.Rotation)]; 
            Offset_TR = InvRot_matrix*Offset_T';
            
            New_OffsetX = mod(Offset_TR(1), PT.XSpacing );
            New_OffsetY = mod(Offset_TR(2), 1/2*PT.YSpacing );
            Rot_matrix = [cos(PT.Rotation),-sin(PT.Rotation);sin(PT.Rotation),cos(PT.Rotation)];
            PT.Offset = (Rot_matrix*[New_OffsetX,New_OffsetY]')';
            
            % shift row estimation 
            NumX = floor(Offset_TR(1)/PT.XSpacing);
            NumY = floor(Offset_TR(2)/(1/2*PT.YSpacing));
            
            if mod(NumX,2) ~= 0 && mod(NumY,2) == 0 || mod(NumX,2) == 0 && mod(NumY,2) ~= 0
            
                PT.FirstPosShiftRow = 1;
                
            else
                PT.FirstPosShiftRow = 2;
            
            end
                  
        end


        function GridCoords = GridBuild(PT)

            % assert number of lenslet in each row is equal, then we can apply same number of lenslet for all rows with
            % confidence
            assert(mod((PT.IM_size(2)-PT.Offset(2)),PT.YSpacing)>PT.YSpacing/2)
            fprintf('Compliant MLA, number of lenslet in each row is equal\n');
            

            Conv2RS_size = PT.IM_size;
            Rot_matrix = [cos(PT.Rotation),-sin(PT.Rotation);sin(PT.Rotation),cos(PT.Rotation)];

            YMax = floor((Conv2RS_size(2)-PT.Offset(2))/(cos(PT.Rotation)*PT.YSpacing))+1 ;
            XMax = floor((Conv2RS_size(1)-PT.Offset(1))/(cos(PT.Rotation)*PT.XSpacing))+1 ;          
            RotTrans_matrix = [Rot_matrix,PT.Offset'];
            [gridX,gridY] = ndgrid((0:XMax-1).*PT.XSpacing, (0:YMax-1).*PT.YSpacing);
            gridY(PT.FirstPosShiftRow:2:end,:) = gridY(PT.FirstPosShiftRow:2:end,:) + 1/2.*PT.YSpacing;
            GridCoords = (RotTrans_matrix*[gridX(:),gridY(:),ones(numel(gridX),1)]')';
            
        end

    end
    
end