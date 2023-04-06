function [filtered_reference, index] = filter_reference(ref, eps)
        %ref is the reference that may include too many solutions, eps is
        %the dimensions of the epsilon box 
        index = 1:size(ref,1);
        for k = 1:length(eps)
            if k~=length(eps)
                [~, idx]=sort(ref(:,k), 'descend');
            else
                [~, idx]=sort(ref(:,k));
            end
            ref = ref(idx,:);
            index = index(idx);
        end
        
        i = 1;
        while i < size(ref,1)
            outOfBox = 0;
            for j = 1:length(eps)
                
                if abs(ref(i, j) - ref(i+1, j)) > eps(j)  %fuori al cubo per soluzione j 
                    outOfBox = 1;
                end
            end
            
            if outOfBox == 0 %soluzione dentro al box, la elimino
                ref(i+1,:) = [];
                index(i+1)   = [];
            else
                i = i+1;
            end
        end
        
        
        filtered_reference = ref;