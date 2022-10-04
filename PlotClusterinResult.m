%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function PlotClusterinResult(X, IDX)
    
    [sizex,sizey] = size(X);
    
    k=max(IDX);

    Colors=lines(k);
    syms = {'*' '.'  'o' 's' '+' 'x' 'd'  '^' 'v' '<' '>' 'p' 'h'};
   % Replicate symbols if number of data points is larger than number of syms
    syms = repmat(syms, 1, ceil((k+1)/numel(syms)));

    %Legends = {};
    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = syms{i+1};
            MarkerSize = 8;
            Color = Colors(i,:);
            %Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = syms{i+1};
            MarkerSize = 6;
            Color = [0 0 0];
            if ~isempty(Xi)
                %Legends{end+1} = 'Noise';
            end
        end
        if ~isempty(Xi)
            if sizey==2
                plot(Xi(:,1),Xi(:,2),Style,'MarkerSize',MarkerSize,'Color',Color);
            elseif sizey == 1
                plot(Xi,zeros(size(Xi)),Style,'MarkerSize',MarkerSize,'Color',Color);
            end
        end
        hold on;
    end
    hold off;
    %axis equal;
    grid on;
    %legend(Legends);
    %legend('Location', 'NorthEastOutside');
end