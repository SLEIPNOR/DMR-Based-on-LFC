function Mesh = OutlierCorrectionEdge(patch_size,binRange,outThres,Mesh,x,y)

Mesh1 = padarray(Mesh, [patch_size, patch_size], 0, 'both');
Mesh1_T = padarray(Mesh, [patch_size, patch_size], 0, 'both');
[zero_y,zero_x] = find(Mesh ~= 0);

% hyper-para setting


for Id= 1:length(zero_y)
    
      
        patch = Mesh1_T(zero_y(Id):zero_y(Id)+patch_size,zero_x(Id):zero_x(Id)+patch_size);
        patch(1/2+size(patch,1)/2,1/2+size(patch,2)/2) = 0;
        patch = patch(:);
        patch = patch(patch~=0);
   
        if ~isempty(patch)
            binEdges = linspace(min(patch),max(patch),(max(patch)-min(patch)+binRange)/binRange);
            estRange = [];
            estStd =[];
        
            if length(binEdges)==1
    
                estStd = mean(patch);
                bias = abs(Mesh1_T(zero_y(Id)+patch_size,zero_x(Id)+patch_size)-estStd);
    
                if bias < outThres 
        
                else
                    
                    Mesh1(zero_y(Id)+patch_size,zero_x(Id)+patch_size) = estStd;
        
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
        
                         if sortedCount(2)>=2 || ismember(Id,[2,x+1]) || ismember(Id,[2,y+1])
    
                            maxBinId_2 = find(histcount == sortedCount(2)); 
        
                            for boundNum = 1:length(maxBinId_2)
        
                                estStd(boundNum+1) = mean(patch(bin==maxBinId_2(boundNum)));
        
                            end
        
                         end          
    
                    end
    
                end
    
                
                
                bias = abs(Mesh1_T(zero_y(Id)+patch_size,zero_x(Id)+patch_size)-estStd);
    
                if any(bias > outThres) 
                    Mesh1(zero_y(Id)+patch_size,zero_x(Id)+patch_size) = mean(estStd); %(bias==min(bias));
        
                else
                    
                    
        
                end
            end
        else
            Mesh1(zero_y(Id)+patch_size,zero_x(Id)+patch_size) = Mesh1(zero_y(Id)+patch_size,zero_x(Id)+patch_size);
        end

    
end
Mesh = Mesh1(1+patch_size:x+patch_size,1+patch_size:y+patch_size);