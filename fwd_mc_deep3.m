function [outcount,outR,outrrs] = fwd_mc_deep3(nzen,nazi,aw,bw,ag,bp,ccgam,cccdf,squad,drad_ang,drad_nph,pcount)
% MC simulation for one IOP set, and one wavelength
% Keep track of depth and angles, but not x-y position
% Runs nph number of photons and outputs no of photons per quad exiting surface

% Initiate outrad
outrad = zeros(nzen,nazi);
outcount = zeros(nzen,nazi);

% Calculate IOPs
at = aw + ag;
bt = bw + bp;
c = at + bt;

% Initiate photon runs
for iang = 1:121 % incoming photon angle
    for iph = 1:drad_nph(iang) % no of photons from iang
       
        % Initialize photon
        outph = 1;
        depth = 0;
        wzen = drad_ang(iang,1);
        wazi = drad_ang(iang,2);
        
        mux = sin(wzen)*cos(wazi);
        muy = sin(wzen)*sin(wazi);
        muz = cos(wzen);
        
        x = (0.55*rand)-0.275;
        y = (0.55*rand)-0.275;
        
        while outph == 1 % photon is within working depth and not absorbed/exited
            
            % Determine path and depth travelled by photon
            path = -log(rand)/c;
            depth = depth + path*muz;
            x = x + path*mux;
            y = y + path*muy;
            
            if depth < 0 % photon exited from water surface
                outph = 2;
                
            elseif (depth > 0.95) || (abs(x) > 0.275) || (abs(y) > 0.275)  % photon hit bottom or side
                if rand < 0.011 % Reflect from the walls with rho = 0.011
                    
                    if depth > 0.95
                        muz = -muz;
                        depth = 0.95;            
                    elseif abs(x) > 0.275
                       muz = -muz;
                       mux = -mux;
                       x = sign(x)*0.275;
                    elseif abs(y) > 0.275
                        muz = -muz;
                        muy = -muy;
                        y = sign(y)*0.275;
                    end
                    
                else
                    outph = 0;
                end
                
            else % still inside working depth
                % Determine if photon is scattered or absorbed
                if rand < (at/c) % absorbed
                    outph = 0;
                    
                else % scattered, outph remains at 1
                    % Determine if scattered by water or particle
                    if rand < (bw/bt) % water
                        q = rand;
                        delzen = acos((2*2^(1/3)*0.835 - 2^(2/3)*(sqrt(0.835^3*(4 + 0.835*(3+0.835)^2*(1-2*q)^2)) + (0.835^2)*(3+0.835)*(-1+2*q))^(2/3))/(2*0.835*(sqrt(0.835^3*(4 + 0.835*(3+0.835)^2*(1-2*q)^2)) + (0.835^2)*(3+0.835)*(-1+2*q))^(1/3)));
                        
                    else % particle
                        [~,ind] = min(abs(cccdf-rand));
                        delzen = ccgam(ind);
                    end
                    
                    % scattering azi
                    delazi = 2*pi*rand;
                    
                    % Calculate new directional cosines after scattering
                    mus = cos(delzen);
                    if abs(muz) < 0.99999
                        muxp = (mux*muz/sqrt(1-muz^2))*(cos(delazi)*sqrt(1-mus^2)) + (-muy/sqrt(1-muz^2))*(sin(delazi)*sqrt(1-mus^2)) + mux*mus;
                        muyp = (muy*muz/sqrt(1-muz^2))*(cos(delazi)*sqrt(1-mus^2)) + (mux/sqrt(1-muz^2))*(sin(delazi)*sqrt(1-mus^2)) + muy*mus;
                        muzp = -sqrt(1-muz^2)*(cos(delazi)*sqrt(1-mus^2)) + muz*mus;
                    else
                        muxp = sign(muz)*cos(delazi)*sqrt(1-mus^2);
                        muyp = sign(muz)*sin(delazi)*sqrt(1-mus^2);
                        muzp = sign(muz)*mus;
                    end
                    
                    % Update directional cosines
                    mux = muxp;
                    muy = muyp;
                    muz = muzp;
                    
                end
                
            end
            
        end
        
        if outph == 2 % water-leaving
            
            % Determine angles (radians) from directional cosines
            outzen = acos(muz);
            
            if muz == -1 % If photon is exiting vertically (exactly 180), assign any azimuth
                muz = -0.9;
            end
            
            if muy >= 0
                outazi = acos(mux/(sqrt(1-muz^2)));
            else
                outazi = 2*pi - acos(mux/(sqrt(1-muz^2)));
            end
            
            
            % Determine which 'quad' the photon came out in (degrees)
            qzen = round(((outzen*180/pi)/10 - 8));
           
            if outazi > pi % Reflect all to one side
                outazi = 2*pi - outazi;
            end
            qazi = round((outazi*180/pi)/15) + 1;

            % Add photon to output storage
            outrad(qzen,qazi) = outrad(qzen,qazi) + 1/abs(muz);
            outcount(qzen,qazi) = outcount(qzen,qazi) + 1;
        end
    end
end

% Sum all at (zen = 180 degrees) into single polar cap
outrad(10,1) = sum(outrad(10,:));
outrad(10,2:end) = 0;

% Calculate R
outR = sum(sum(outcount))/(pcount);

% Calculate rrs
outrrs = outrad./(2*squad*pcount);

end
