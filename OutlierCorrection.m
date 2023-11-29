function Mesh = OutlierCorrection(patch_size,binRange,outThres,Mesh,x,y)

Mesh1 = padarray(Mesh, [patch_size, patch_size], 0, 'both');
Mesh2 = Mesh1;
Mesh3 = Mesh1;
Mesh4 = Mesh1;

% hyper-para setting


for Id_row = 1+patch_size:x+patch_size
    for Id_col = 1+patch_size:y+patch_size %1+1:n+1
        
        patch = Mesh1(Id_row-patch_size:Id_row+patch_size,Id_col-patch_size:Id_col+patch_size);
        patch(1/2+size(patch,1)/2,1/2+size(patch,2)/2) = 0;
        patch = patch(:);
        patch = patch(patch~=0);
        binEdges = linspace(min(patch),max(patch),(max(patch)-min(patch)+binRange)/binRange);
        estRange = [];
        estStd =[];
        if length(binEdges)==1

            estStd = mean(patch);
            bias = abs(Mesh1(Id_row,Id_col)-estStd);

            if bias < outThres 
    
            else
                
                Mesh1(Id_row,Id_col) = estStd;
    
            end

        else

            [histcount,binEdges,bin]=histcounts(patch,binEdges);

            maxBinId = find(histcount == max(histcount));
            if length(maxBinId) > 1

                for boundNum = 1:length(maxBinId)
        
                      estRange = [estRange, patch(bin==maxBinId(boundNum))];  
    
                end

                estStd = mean(estRange,1);

            else
                estRange = patch(bin==maxBinId);
                estStd(1) = mean(estRange);

                if length(histcount)>1

                     sortedCount = sort(unique(histcount),'descend');
    
                     if sortedCount(2)>=2 || ismember(Id_row,[2,x+1]) || ismember(Id_col,[2,y+1])

                        maxBinId_2 = find(histcount == sortedCount(2)); 
    
                        for boundNum = 1:length(maxBinId_2)
    
                            estStd(boundNum+1) = mean(patch(bin==maxBinId_2(boundNum)));
    
                        end
    
                     end          

                end

            end

            
            
            bias = abs(Mesh1(Id_row,Id_col)-estStd);

            if any(bias < outThres) 
    
            else
                
                Mesh1(Id_row,Id_col) = mean(estStd); %(bias==min(bias));
    
            end
        end

    end
end

Mesh1 = Mesh1(1+patch_size:x+patch_size,1+patch_size:y+patch_size);

%**************************************************************************

for Id_row = x+patch_size:-1:1+patch_size
    for Id_col = 1+patch_size:y+patch_size %1+1:n+1

        patch = Mesh2(Id_row-patch_size:Id_row+patch_size,Id_col-patch_size:Id_col+patch_size);
        patch(1/2+size(patch,1)/2,1/2+size(patch,2)/2) = 0;
        patch = patch(:);
        patch = patch(patch~=0);
        binEdges = linspace(min(patch),max(patch),(max(patch)-min(patch)+binRange)/binRange);
        estRange = [];
        estStd =[];
        if length(binEdges)==1

            estStd = mean(patch);
            bias = abs(Mesh2(Id_row,Id_col)-estStd);

            if bias < outThres 
    
            else
                
                Mesh2(Id_row,Id_col) = estStd;
    
            end

        else

            [histcount,binEdges,bin]=histcounts(patch,binEdges);

            maxBinId = find(histcount == max(histcount));
            if length(maxBinId) > 1

                for boundNum = 1:length(maxBinId)
        
                      estRange = [estRange, patch(bin==maxBinId(boundNum))];  
    
                end

                estStd = mean(estRange,1);

            else
                estRange = patch(bin==maxBinId);
                estStd(1) = mean(estRange);

                if length(histcount)>1

                     sortedCount = sort(unique(histcount),'descend');
    
                     if sortedCount(2)>=2 || ismember(Id_row,[2,x+1]) || ismember(Id_col,[2,y+1])

                        maxBinId_2 = find(histcount == sortedCount(2)); 
    
                        for boundNum = 1:length(maxBinId_2)
    
                            estStd(boundNum+1) = mean(patch(bin==maxBinId_2(boundNum)));
    
                        end
    
                     end          

                end

            end

            
            
            bias = abs(Mesh2(Id_row,Id_col)-estStd);

            if any(bias < outThres) 
    
            else
                
                Mesh2(Id_row,Id_col) = mean(estStd); %(bias==min(bias));
    
            end
        end

    end
end

Mesh2 = Mesh2(1+patch_size:x+patch_size,1+patch_size:y+patch_size);
%**************************************************************************

for Id_row = 1+patch_size:x+patch_size
    for Id_col = y+patch_size:-1:1+patch_size %1+1:n+1

        patch = Mesh3(Id_row-patch_size:Id_row+patch_size,Id_col-patch_size:Id_col+patch_size);
        patch(1/2+size(patch,1)/2,1/2+size(patch,2)/2) = 0;
        patch = patch(:);
        patch = patch(patch~=0);
        binEdges = linspace(min(patch),max(patch),(max(patch)-min(patch)+binRange)/binRange);
        estRange = [];
        estStd =[];
        if length(binEdges)==1

            estStd = mean(patch);
            bias = abs(Mesh3(Id_row,Id_col)-estStd);

            if bias < outThres 
    
            else
                
                Mesh3(Id_row,Id_col) = estStd;
    
            end

        else

            [histcount,binEdges,bin]=histcounts(patch,binEdges);

            maxBinId = find(histcount == max(histcount));
            if length(maxBinId) > 1

                for boundNum = 1:length(maxBinId)
        
                      estRange = [estRange, patch(bin==maxBinId(boundNum))];  
    
                end

                estStd = mean(estRange,1);

            else
                estRange = patch(bin==maxBinId);
                estStd(1) = mean(estRange);

                if length(histcount)>1

                     sortedCount = sort(unique(histcount),'descend');
    
                     if sortedCount(2)>=2 || ismember(Id_row,[2,x+1]) || ismember(Id_col,[2,y+1])

                        maxBinId_2 = find(histcount == sortedCount(2)); 
    
                        for boundNum = 1:length(maxBinId_2)
    
                            estStd(boundNum+1) = mean(patch(bin==maxBinId_2(boundNum)));
    
                        end
    
                     end          

                end

            end

            
            
            bias = abs(Mesh3(Id_row,Id_col)-estStd);

            if any(bias < outThres) 
    
            else
                
                Mesh3(Id_row,Id_col) = mean(estStd); %(bias==min(bias));
    
            end
        end

    end
end

Mesh3 = Mesh3(1+patch_size:x+patch_size,1+patch_size:y+patch_size);

%**************************************************************************

for Id_row = x+patch_size:-1:1+patch_size
    for Id_col = y+patch_size:-1:1+patch_size %1+1:n+1

        patch = Mesh4(Id_row-patch_size:Id_row+patch_size,Id_col-patch_size:Id_col+patch_size);
        patch(1/2+size(patch,1)/2,1/2+size(patch,2)/2) = 0;
        patch = patch(:);
        patch = patch(patch~=0);
        binEdges = linspace(min(patch),max(patch),(max(patch)-min(patch)+binRange)/binRange);
        estRange = [];
        estStd =[];
        if length(binEdges)==1

            estStd = mean(patch);
            bias = abs(Mesh4(Id_row,Id_col)-estStd);

            if bias < outThres 
    
            else
                
                Mesh4(Id_row,Id_col) = estStd;
    
            end

        else

            [histcount,binEdges,bin]=histcounts(patch,binEdges);

            maxBinId = find(histcount == max(histcount));
            if length(maxBinId) > 1

                for boundNum = 1:length(maxBinId)
        
                      estRange = [estRange, patch(bin==maxBinId(boundNum))];  
    
                end

                estStd = mean(estRange,1);

            else
                estRange = patch(bin==maxBinId);
                estStd(1) = mean(estRange);

                if length(histcount)>1

                     sortedCount = sort(unique(histcount),'descend');
    
                     if sortedCount(2)>=2 || ismember(Id_row,[2,x+1]) || ismember(Id_col,[2,y+1])

                        maxBinId_2 = find(histcount == sortedCount(2)); 
    
                        for boundNum = 1:length(maxBinId_2)
    
                            estStd(boundNum+1) = mean(patch(bin==maxBinId_2(boundNum)));
    
                        end
    
                     end          

                end

            end

            
            
            bias = abs(Mesh4(Id_row,Id_col)-estStd);

            if any(bias < outThres) 
    
            else
                
                Mesh4(Id_row,Id_col) = mean(estStd); %(bias==min(bias));
    
            end
        end

    end
end

Mesh4 = Mesh4(1+patch_size:x+patch_size,1+patch_size:y+patch_size);

%**************************************************************************
Mesh = round((Mesh1+Mesh2+Mesh3+Mesh4)/4);
