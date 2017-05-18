function view3D(W,S, edges)

    F = size(W,1)/2;
    
    for i=1:F 

        plot3(S(3*i-2, :), S(3*i, :), S(3*i-1, :), 'k.'); 
        
        drawLines(S(3*i-2:3*i, :), 'k-',edges);
       
       % hold off;
       % title(sprintf('Frame: %03d',i))
        title('3D Human Poses')
        set(gca,'xaxislocation','top','yaxislocation','left','zdir','reverse','xdir','reverse')
        view(180, 24); 
        
        grid on
      
        axis equal
        set(gca,'xtick',[],'ytick',[],'ztick',[]);
%        
    end  
        
end







function drawLines(S,lineStyle,edges)
    
    for b = 1:size(edges,1)
        i = edges(b,1); j = edges(b,2);
        plot3([S(1,i), S(1,j)], [S(3,i), S(3, j)], [S(2,i), S(2, j)], lineStyle, 'linewidth', 1.5);            
    end        
end