function TTLstartsTimes=GetTTLtimes(TTLData,ExtractedResults)

dTTL = diff(TTLData);                                                             % Time derivation of the TTL signal to catch TTL rises

dTTL(dTTL<max(dTTL)/50) = 0;                                                    % Set to 0 everything below a half of the max
TTLstarts = find(dTTL);                                                         % Find the position of every remaining non zero values (strong rises)

distStarts = diff(TTLstarts);                               % find distance between putative rises (should be around 200 because TTL is 100Hz so 1 every 10ms=200pts at 20kHz)
TTLstarts(distStarts < 0.5*median(distStarts)) = [];        % remove points that are too close (<medianDist/2 ?) (Happen when rises is on 3 points --> 2 points kept)

distStarts = diff(TTLstarts);                               % Update distances (it should be a while loop)
TTLBreaks = find(distStarts > 3*median(distStarts));        % Find gap between TTLs, i.e when TTLs are spaced by more than 3 times the median
TTLBreaks=[TTLBreaks(1) TTLBreaks(end)];                    % Keep the first and the last value (just in case there is a gap in the middle, hopefully there is no gap before the start and after the end
TTLBreaksTimes = [(TTLstarts(TTLBreaks(1))+TTLstarts(TTLBreaks(1)+1))/2 (TTLstarts(TTLBreaks(2))+TTLstarts(TTLBreaks(2)+1))/2]; % Timings of the breaks

TTLstarts((TTLBreaks(2)+1):end) = [];                       % remove TTL after second break (so after the end of the experiment)
TTLstarts(1:TTLBreaks(1)) = [];                             % remove TTL before firs break (so before the beginning of the experiment)

disp(['4 firsts TTL were removed ; ' num2str(numel(TTLstarts)-4-numel(ExtractedResults.C_raw(1,:))) ' were removed at the end'])
TTLstarts=TTLstarts(5 : 4+numel(ExtractedResults.C_raw(1,:)));       % remove the 4 first TTL (images are discarded) and any surnumerary TTLs at the end

TTLstartsTimes=TTLstarts/20000;