function view2D(W, edges)

    F = size(W,1)/2;

    for i=1:F 
              
        plot(W(2*i-1, :), W(2*i, :), 'k.'); 
          
        drawLines2(W(2*i-1:2*i, :), 'k-',edges);        
        title('2D Skeletos')
        grid on
        set(gca,'xaxislocation','top','yaxislocation','left','ydir','reverse')
        axis equal
    end   
     
end

function drawLines2(W,lineStyle,edges)
    
    for b = 1:size(edges,1)
        i = edges(b,1); j = edges(b,2);
        plot([W(1,i), W(1,j)], [W(2,i), W(2, j)], lineStyle);            
    end        
end


