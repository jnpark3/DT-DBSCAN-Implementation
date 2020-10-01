% MATLAB code for triangulation visualization
% Author: Jian Park (21PARKJ@sagehillschool.org)

visnow( 'trig.3.out', -1, -1 );

function circle(x,y,r)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,'.-');
end

function visnow( dsource, cidx, r )
    fid = fopen( dsource );
    line = fgets( fid );
    N = str2double(line);
    X = linspace(0,0,N);
    Y = linspace(0,0,N);

    for i = 1:N
        line = fgets( fid );
        X(i) = str2double(line);

        line = fgets( fid );
        Y(i) = str2double(line);
    end
    
    xMin = min(X);
    xMax = max(X);
    yMin = min(Y);
    yMax = max(Y);
    gridWidth = max( (yMax-yMin), (xMax-xMin) );
    k1 = ( uint32( log10(N)/log10(4) ) );
    K = 4 ^ k1;
    gridnum = sqrt( double(K) );
    gridsize = gridWidth / gridnum;

    for x = xMin:gridsize:max(xMax,yMax)
        GX = linspace(0,0,2);
        GY = linspace(0,0,2);
        GX(1) = x;
        GY(1) = yMin;
        GX(2) = x;
        GY(2) = yMax;
        plot(GX,GY,'Color',[0.05,0.05,0.05]);
        hold on
    end

    for y = yMin:gridsize:max(xMax,yMax)
        GX = linspace(0,0,2);
        GY = linspace(0,0,2);
        GX(1) = xMin;
        GY(1) = y;
        GX(2) = xMax;
        GY(2) = y;
        plot(GX,GY,'Color',[0.05,0.05,0.05]);
        hold on
    end

    while (1)
        line = fgets( fid );
        if ( line == - 1 )
            break;
        end
        vIdx = strsplit( line );
        cnt = numel(vIdx)-1;
        
        sIdx = str2double( vIdx(1) );
        if ( sIdx < 0 )
            sIdx = (-1) * sIdx;
            % text(X(sIdx),Y(sIdx),'.','Color','red');
            circle(X(sIdx),Y(sIdx),5);
        end
        sIdx = sIdx + 1;

        VX = linspace(0,0,cnt);
        VY = linspace(0,0,cnt);
        vCnt = 1;
        VX(vCnt) = X(sIdx);
        VY(vCnt) = Y(sIdx);

%        text(X(sIdx),Y(sIdx),sprintf("%d",sIdx-1),'Color',color);

        for i = 2:cnt
            % color = 'white';
            nIdx = str2double( vIdx(i) );
            if ( nIdx < 0 )
                nIdx = (-1) * nIdx;
                % color = 'red';
                % text(X(nIdx),Y(nIdx),'.','Color','red');
            end
            nIdx = nIdx + 1;
            
            VX(vCnt+1) = X(nIdx);
            VY(vCnt+1) = Y(nIdx);
            VX(vCnt+2) = X(sIdx);
            VY(vCnt+2) = Y(sIdx);
            vCnt = vCnt + 2;

            % text(X(nIdx),Y(nIdx),sprintf("%d",nIdx-1),'Color',color);
        end
        plot(VX,VY,'.-', 'Color', [0 102 204]/256 );
        hold on;
    end
    
    if ( cidx >= 0 )
        circle( X( cidx + 1 ), Y( cidx + 1 ), r );
    end
    
    outfile = strcat( dsource, '.png' );
    saveas(gcf, outfile );
end
