function Hs = ChannelUpdate(ChObjects,numUEs,Nsubray,BSobj,UEobj)
Mt = BSobj.M; Nt = BSobj.N;
Mr = UEobj.M; Nr = UEobj.N;
% c  = 299792458; 
f  = 28e9; %lamda = c/f;
Hs = zeros(numUEs,Mt*Nt,Mr*Nr);
%Orientation = BSobj.Orientation;
for  UEnow  = 1:numUEs
    Object  = ChObjects{UEnow};
    N = length(Object);
    H = zeros(Mt*Nt,Mr*Nr);
    for raynow = 1:N
        %==========some randomness here ===================================
        if ChObjects{UEnow}{raynow}.LOS==0
            DoDs_phy   = ChObjects{UEnow}{raynow}.DoDs_phy   + unifrnd(0,0.01/180*pi);  % update of sub-ray angles
            DoDs_theta = ChObjects{UEnow}{raynow}.DoDs_theta + unifrnd(0,0.01/180*pi);
            DoAs_phy   = ChObjects{UEnow}{raynow}.DoAs_phy   + unifrnd(0,0.01/180*pi);
            DoAs_theta = ChObjects{UEnow}{raynow}.DoAs_theta + unifrnd(0,0.01/180*pi); 
            
            %meandelay = 10*1e-9; % the mean delay in a cluster associated with a cluster
            %subrayDelays = exprnd(meandelay,1,Nsubray);                                % update of sub-ray delays
            NNsubray = Nsubray;
        else
            DoDs_phy = 0;
            DoDs_theta = 0;
            DoAs_phy = 0;
            DoAs_theta = 0;
            %subrayDelays = 0;
            %subrayPowers = 1;
            NNsubray = 1; % LOS path,one sub-ray only
        end
        
        %ChObjects{UEnow}{raynow}.DoDs_phy   = DoDs_phy;
        %ChObjects{UEnow}{raynow}.DoDs_theta = DoDs_theta;
        %ChObjects{UEnow}{raynow}.DoAs_phy   = DoAs_phy;
        %ChObjects{UEnow}{raynow}.DoAs_theta = DoAs_theta;
        
        %ChObjects{UEnow}{raynow}.subrayDelays = subrayDelays;  % update of sub-ray delays
        
        % keep the subray powers
        subrayPowers = ChObjects{UEnow}{raynow}.subrayPowers;

        % Cluster
        % see: [1] 28 GHz Channel Modeling Using 3D Ray-tracing in Uran Environments
        for subray = 1:NNsubray
            
            %aBS   = arrayresponse(Object{raynow}.DoD + DoDs_theta(subray) - Orientation, DoDs_phy(subray), BSobj);
            aBS   = arrayresponse(Object{raynow}.DoD + DoDs_theta(subray), DoDs_phy(subray), BSobj);
            aUE   = arrayresponse(Object{raynow}.DoA + DoAs_theta(subray), DoAs_phy(subray), UEobj);
            
            %=====Patched antenna==========================================
            %ElementGainBS = ElementPattern( Object{raynow}.DoD + DoDs_theta(subray) - Orientation, DoDs_phy(subray), BSobj);
            ElementGainBS = ElementPattern(Object{raynow}.DoD + DoDs_theta(subray), DoDs_phy(subray), BSobj);
            ElementGainUE = ElementPattern(Object{raynow}.DoA + DoAs_theta(subray), DoAs_phy(subray), UEobj);
            %3-D Statistical Channel Model for Millimeter-Wave Outdoor Mobile Broadband Communications
            subrayDelays = ChObjects{UEnow}{raynow}.subrayDelays;
            subpathdelay = Object{raynow}.delay + subrayDelays(subray);
            subpathphase = subpathdelay * f * 2 * pi;
            alpha = exp(2*1j*pi*rand)*ElementGainBS*ElementGainUE*(Object{raynow}.amplitude*sqrt(subrayPowers(subray)))*exp(1i*subpathphase);  % random phase of the ray arrived 
            
            
            H = H + alpha*aBS*aUE';         %////////////
        end
    end
    Hs(UEnow,:,:) = H;
end
end