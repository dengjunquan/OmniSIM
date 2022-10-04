function [Hs, ChObjects]= SalehValenzuela(path_arrays,walls,eta_rs,sigma_theta,sigma_phy,max_theta,max_phy,Nsubray,BSobj,UEobj,Orientation)
% for each BS-to-UE links
% we have N rays
% we have N DoDs and N DoAs
% we have N ray delays
% we have N ray amplitudes dependent on the random
% reflection coefficients(related to incident angles)

c  = 299792458; f  = 28e9; %lamda = c/f;
numUEs = length(path_arrays);
ChObjects  = cell(numUEs,1);
gg = 1/sqrt(10^(61.4/10));

zeroface = zeros(1,4);
nface = zeros(1,4);

for UEnow  = 1:numUEs
    path_array = path_arrays{UEnow};
    N = length(path_array); % number of rays
    Object = cell(N,1);
    if isempty(path_array)
        ChObjects{UEnow} = Object;
        continue;
    end
    
    anypath = path_array{1};
    UExy = anypath(size(anypath,1),1:2); 
    BSxy = anypath(1,1:2);  
    
    for raynow = 1:N % for each ray (cluster)
        
        pathnow = path_array{raynow};
        % 多径类型 纯反射或包含衍射
        pathtype  = pathnow(1,end); % 1 reflection, 2 diffraction
        
        % DoD at BS
        StartXY = BSxy;  % equals pathnow(1,1:2)
        EndedXY = pathnow(2,1:2);
        Object{raynow}.DoD = atan2(EndedXY(2)-StartXY(2),EndedXY(1)-StartXY(1)) - Orientation; % atan2(y,x)  -pi to pi
        
        %disp(['UE = ',num2str(UEnow),',',num2str(Object{raynow}.DoD*180/pi)]);
        
        % DoA at UE
        EndedXY = pathnow( size(pathnow,1) - 1,1:2 );
        StartXY = UExy;   % equals pathnow(end,1:2)
        Object{raynow}.DoA = atan2(EndedXY(2)-StartXY(2),EndedXY(1)-StartXY(1)); % atan2(y,x)  -pi to pi
        
        if  size(pathnow,1)==2  % This is a LOS ray
            Object{raynow}.LOS = 1;
            Object{raynow}.order = 0;
            distance = sqrt((BSxy(1)-UExy(1))^2 + (BSxy(2)-UExy(2))^2);
            d_dash = distance; 
            Object{raynow}.delay = distance/c;
            Object{raynow}.amplitude = gg*(1/distance);
            
        else  % a multipath ray
            Object{raynow}.LOS = 0;
            Object{raynow}.order = size(pathnow,1)-2;
            power_ratio = 1;
            
            StartXY = pathnow(1,1:2);
            EndedXY = pathnow(2,1:2);     
            % a path is represented by 
            % BSxy <- point1 <-point2 <-...<-pointN-1 <- UExy 
            % with Nsegments
            
            distance = sqrt((EndedXY(1)-StartXY(1))^2+(EndedXY(2)-StartXY(2))^2);
            d_dash = distance; 
            
            for Pointnow = 2 : size(pathnow,1)-1
                
                StartXY = pathnow(Pointnow,    1:2); 
                EndedXY = pathnow(Pointnow + 1,1:2); 
                d = sqrt((EndedXY(1)-StartXY(1))^2+(EndedXY(2)-StartXY(2))^2);
                distance = distance + d;
                
                if Pointnow ==2 && pathtype == 2 % diffraction
                    
                    Wedge = StartXY;                       % diffraction point
                    zerofaceID = pathnow(Pointnow,3);      % ID for the 0-face wall
                    nfaceID    = pathnow(Pointnow,4);      % ID for the n-face wall
                    
                    zeroface(1) = walls{zerofaceID,1}.p1(1);
                    zeroface(2) = walls{zerofaceID,1}.p1(2);
                    zeroface(3) = walls{zerofaceID,1}.p2(1);
                    zeroface(4) = walls{zerofaceID,1}.p2(2);
                    
                    nface(1) = walls{nfaceID,1}.p1(1);
                    nface(2) = walls{nfaceID,1}.p1(2);
                    nface(3) = walls{nfaceID,1}.p2(1);
                    nface(4) = walls{nfaceID,1}.p2(2);
                    
                    dX = Wedge(1)-BSxy(1); dY = Wedge(2)-BSxy(2);
                    phi_in = atan2(dY,dX);
                    
                    % 判断衍射体的边界朝向
                    if nface(1) == Wedge(1)
                        dXn = nface(3)-nface(1); dYn = nface(4)-nface(2);
                    else
                        dXn = nface(1)-nface(3); dYn = nface(2)-nface(4);
                    end
                    
                    if zeroface(3) == Wedge(1)
                        dX0 = zeroface(3)-zeroface(1); dY0 = zeroface(4)-zeroface(2);
                    else
                        dX0 = zeroface(1)-zeroface(3); dY0 = zeroface(2)-zeroface(4);
                    end
                    
                    v_out = [EndedXY(1) - StartXY(1), EndedXY(2)- StartXY(2)];
                    phi_out = atan2(v_out(2),v_out(1));
                    
                    phi_nf = atan2(dYn,dXn);
                    phi_0f = atan2(dY0,dX0);
                    
                    phi_in  = phi_in  + (phi_in <0)*2*pi; % 入射波方向
                    phi_out = phi_out + (phi_out<0)*2*pi; % 出射波方向
                    
                    phi = abs(phi_in-phi_out);
                    
                    if phi < pi
                        phi = pi + phi;
                    end
                    
                    phi_nf = phi_nf + (phi_nf<0)*2*pi; % n-face 方向
                    phi_0f = phi_0f + (phi_0f<0)*2*pi; % zero-face 方向
                    
                    phi_dash = abs(phi_in-phi_0f);     % 入射波与zero-face夹角
                    
                    if phi_dash>pi
                        phi_dash = 2*pi - phi_dash;
                    end
                    
                    % A New Heuristic UTD Diffraction Coefficient for Nonperfectly Conducting Wedges
                    eta_0 = eta_rs(zerofaceID);
                    eta_n = eta_rs(nfaceID);
                    phi_0n = abs(phi_nf-phi_0f);
                    if phi_0n>pi
                        phi_0n = 2*pi - phi_0n;
                    end
                    n = (2*pi - phi_0n)/pi;
                    
                    ToPlot = 0;
                    if ToPlot == 1
                        figure;hold on;
                        scatter(Wedge(1),Wedge(2),'k*');
                        scatter(BSxy(1),BSxy(2),'b^');
                        scatter(StartXY(1),StartXY(2),'rs');
                        scatter(EndedXY(1),EndedXY(2),'rd');
                        scatter(UExy(1),UExy(2),'bO');
                        plot(zeroface(1:2:3),zeroface(2:2:4),'b-');
                        plot(nface(1:2:3),nface(2:2:4),'b-');
                        axis equal;
                    end
    
                    diffcoeff = DiffractionCoeff(n,d_dash,d,phi_dash,phi,eta_0,eta_n,f);
                    power_ratio  = power_ratio*diffcoeff^2;
                    
                else
                    wallid  = pathnow(Pointnow,3);       % ID for the wall
                    v_out   = [EndedXY(1) - StartXY(1), EndedXY(2)- StartXY(2),0];
                    v_wall  = [walls{wallid}.p1(1) - walls{wallid}.p2(1), walls{wallid}.p1(2)-walls{wallid}.p2(2),0];
                    theta   = atan2(norm(cross(v_out,v_wall)),dot(v_out,v_wall));
                    if theta > pi/2
                        theta = pi - theta;   % angle between the reflected ray and the wall face
                    end
                    theta_i = pi/2 - theta;   % the incident angele
                    
                    coeff        = ReflectionCoeff(theta_i,eta_rs(wallid));
                    power_ratio  = power_ratio*coeff^2;
                end
            end
            
            Object{raynow}.distanceAfterInteraction = distance - d_dash;
            Object{raynow}.delay = distance/c;
            Object{raynow}.amplitude = gg*(1/distance)*sqrt(power_ratio);
        end
    end
    ChObjects{UEnow} = Object;
end

Mt = BSobj.M; Nt = BSobj.N;
Mr = UEobj.M; Nr = UEobj.N;

Hs = zeros(numUEs,Mt*Nt,Mr*Nr);

meandelaymin = 0.5/c;  %10*1e-9;  % the mean delay in a cluster associated with a cluster
            
for  UEnow  = 1:numUEs
    Object  = ChObjects{UEnow};
    N = length(Object);
    H = zeros(Mt*Nt,Mr*Nr);
    for raynow = 1:N
        %==========some randomness here ===================================
        if Object{raynow}.LOS==0
            pathorder  = ChObjects{UEnow}{raynow}.order;
            %distanceAfterInteraction = Object{raynow}.distanceAfterInteraction;
            DoDs_phy   = laprnd(Nsubray,0,sigma_phy*pathorder,max_phy*pathorder);
            DoDs_theta = laprnd(Nsubray,0,sigma_theta*pathorder,max_theta*pathorder);
            DoAs_phy   = laprnd(Nsubray,0,sigma_phy*pathorder,max_phy*pathorder);
            DoAs_theta = laprnd(Nsubray,0,sigma_theta*pathorder,max_theta*pathorder);
            
            %===========================================================
            % Delays for cluster subpaths are exponentially distributed
            meandelay = meandelaymin*Object{raynow}.order;
            subrayDelays = exprnd(meandelay,1,Nsubray);
            %subrayDelays = sort(subrayDelays);
            
            % 3-D Statistical Channel Model for Millimeter-Wave Outdoor Mobile Broadband Communications
            subrayPowers = exp(-subrayDelays./meandelay).*10.^(randn(1,Nsubray)./10);
            subrayPowers = subrayPowers./sum(subrayPowers);
            
            subrayPowers = sort(subrayPowers,'descend');
            
            % Mapping between powers, delays, angles 
            
            %[~,IndexDoD_theta] = sort(abs(DoDs_theta));
            %[~,IndexDoA_theta] = sort(abs(DoAs_theta));
            %DoDs_theta = DoDs_theta(IndexDoD_theta);
            %DoAs_theta = DoAs_theta(IndexDoA_theta);
            
            %[~,IndexDoD_phy] = sort(abs(DoDs_phy));
            %[~,IndexDoA_phy] = sort(abs(DoAs_phy));
            %DoDs_phy = DoDs_theta(IndexDoD_phy);
            %DoAs_phy = DoAs_theta(IndexDoA_phy);
            
            % plot(subrayDelays,subrayPowers,'O');
            %===========================================================
            NNsubray = Nsubray;
        else
            DoDs_phy = 0;
            DoDs_theta = 0;
            DoAs_phy = 0;
            DoAs_theta = 0;
            subrayDelays = 0;
            subrayPowers = 1;
            NNsubray = 1; % LOS path,one sub-ray only
        end
        
        ChObjects{UEnow}{raynow}.DoDs_phy   = DoDs_phy;
        ChObjects{UEnow}{raynow}.DoDs_theta = DoDs_theta;
        ChObjects{UEnow}{raynow}.DoAs_phy   = DoAs_phy;
        ChObjects{UEnow}{raynow}.DoAs_theta = DoAs_theta;
        ChObjects{UEnow}{raynow}.subrayDelays = subrayDelays;
        ChObjects{UEnow}{raynow}.subrayPowers = subrayPowers;
        
        ElementGainBS = ElementPattern( Object{raynow}.DoD, 0, BSobj);
        ElementGainUE = ElementPattern( Object{raynow}.DoA, 0, UEobj);
        ChObjects{UEnow}{raynow}.ClusterPower = (ElementGainBS*ElementGainUE*ChObjects{UEnow}{raynow}.amplitude)^2;
        
        % Cluster
        % see: [x] 28 GHz Channel Modeling Using 3D Ray-tracing in Uran Environments
        for subray = 1:NNsubray
            
            aBS   = arrayresponse(Object{raynow}.DoD + DoDs_theta(subray), DoDs_phy(subray), BSobj);
            aUE   = arrayresponse(Object{raynow}.DoA + DoAs_theta(subray), DoAs_phy(subray), UEobj);
            
            % aBS   = arrayresponse(Object{raynow}.DoD, pi/2, Mt,Nt,BSobj.dx,BSobj.dz); % no sub-ray
            % aUE   = arrayresponse(Object{raynow}.DoA, pi/2, Mr,Nr,UEobj.dx,UEobj.dx);
            % the cluster power is equally distributed to the subpath power
            
            %=====Patched antenna==========================================
            ElementGainBS = ElementPattern( Object{raynow}.DoD + DoDs_theta(subray), DoDs_phy(subray), BSobj);
            ElementGainUE = ElementPattern( Object{raynow}.DoA + DoAs_theta(subray), DoAs_phy(subray), UEobj);
            %3-D Statistical Channel Model for Millimeter-Wave Outdoor Mobile Broadband Communications
            subpathdelay = Object{raynow}.delay + subrayDelays(subray);
            subpathphase = subpathdelay * f * 2 * pi;
            alpha = exp(2*1j*pi*rand)*ElementGainBS*ElementGainUE*(Object{raynow}.amplitude*sqrt(subrayPowers(subray)))*exp(1i*subpathphase);  % random phase of the ray arrived
            
            %alpha = (Object{raynow}.amplitude/sqrt(Nsubray))*exp(1i* Object{raynow}.delay * f * 2*pi );
            
            H = H + alpha*aBS*aUE';   %////////////
        end
    end
    Hs(UEnow,:,:) = H;
end
end