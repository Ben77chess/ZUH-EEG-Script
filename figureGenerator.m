% This file is derived from the main script ("test_Conversion_v2.m", currently)
% It generates a number of figures for reader appraisal so that ultimately
% Bland-Altman plots can be created.
% Benjamin Weinberg - Fall 2018

%nicoletFiles = dir('**/*.e');
figurecounter = 2;
%filesRead = strings(1, 16800); %[""] The magic number here is the total number of files being scanned (preallocation hopefully for speed)
%filesError = [""];
%noEyesClosed = [""];
%fileCounter = 0;
%oldFileCounter = 0;
o1pksarr = {}; %Deprecated
o2pksarr = {}; %Deprecated
lowCutoffStatic = 6;
lowCutoff = 6;
highCutoff = 13;
bandpassLowCutoff = .5;
bandpassHighCutoff = 70;
childLowCutoff = .5;
childAge = 8;
secondsToSample = 5; %How many seconds of data do we want to pull?
desiredChannels = ["O1", "O2"];
secondOffset = 1; %Wating this long after eyes closed to collect data
componentGraphs = false; % Do we want to plot each eyes closed event we analyze + corresponding signal?
patientInfo = true; %If we have the info or not (Needs this in order to function in a meaningful capacity)
%dataMap = containers.Map('KeyType','char', 'ValueType','any');
scrambledID = 0;
oldDB = false; %Deprecated
fractionToEval = 1;

ylims = [-50 50]
manualFiles = "";
autoEstimates = "";
while(length(manualFiles) < 11)
    lowCutoff = lowCutoffStatic;
    rndmFName = filesRead(randi(length(filesRead)));
    seg_counter = 1;
    if(~any(contains(manualFiles, rndmFName)) && ~any(contains(noEyesClosed, rndmFName))) %If we haven't sampled it randomly yet and it has EC events
        manualFiles = [manualFiles rndmFName];
        OBJ = NicoletFile(rndmFName); %This part is all going to be duplicated \
        if(patientInfo)
            OBJ.patientInfo.altID(regexp(OBJ.patientInfo.altID,'[-]'))=[];
            birth_dmy = [str2double(OBJ.patientInfo.altID(1:2)), str2double(OBJ.patientInfo.altID(3:4)), str2double(OBJ.patientInfo.altID(5:6))+1900];
            EEG_dmy = [OBJ.segments(1).startDate(3), OBJ.segments(1).startDate(2), OBJ.segments(1).startDate(1)];
            %age = EEG_dmy(3) - birth_dmy(3);
            if (str2double(OBJ.patientInfo.altID(7)) > 4) %Born in the 21st century
                birth_dmy(3) = str2double(OBJ.patientInfo.altID(5:6))+2000;
            end
            %    age = mod(age, 100); %should account for people > 100
            %end
            %if( (EEG_dmy(2) - birth_dmy(2)) > 0 || ( (EEG_dmy(2) - birth_dmy(2)) == 0 && ( (EEG_dmy(1) - birth_dmy(1)) >= 0) ) )
            %    age = age+1;
            %end
            t1 = datetime(birth_dmy(3), birth_dmy(2), birth_dmy(1));
            t2 = datetime(EEG_dmy(3), EEG_dmy(2), EEG_dmy(1));
            age = years(t2-t1);
            sex = mod(str2double(OBJ.patientInfo.altID(10)),2); %1 = M, 0 = F
            if(isnan(age) || isnan(sex))
                ME = MException('Age or Sex is NaN in patient ' + [file.folder '\' file.name]);
                throw(ME)
            end
            scrambledID = base64toHex(doHMAC_SHA1(OBJ.patientInfo.altID, "Changed")); %Change the key in here if the code is ever published
        end
        if(age < childAge)
            lowCutoff = childLowCutoff;
        end
        sampleRate = getSampleRate(OBJ,1,1); %Sample rate of 1st sample of 1st section
        nChans = length(OBJ.segments(1).chName);
        desiredChannelIndices = [];
        for x = 1:length(desiredChannels) %Gets the indices of the channels we want to plot, strcmp returns a logical array with a true in the positon of the channel name, Find retrieves the index of that true
            if find(strcmp(OBJ.segments(1).chName, desiredChannels(x))) ~= 0 %Only checks first segment rn
                desiredChannelIndices = [desiredChannelIndices find(contains(OBJ.segments(1).chName, desiredChannels(x)))];
            end
        end
        
        %Get the indices of all eyes closed events
        eyesClosedEventIndices = [];
        for x = 1:length(OBJ.eventMarkers)
            if OBJ.eventMarkers(x).IDStr == "Eyes Closed"
                eyesClosedEventIndices = [eyesClosedEventIndices x];
            end
        end
        
        for x = randi(length(eyesClosedEventIndices))
            try
                tindx = getTimeIndexOfEvent(OBJ, eyesClosedEventIndices(x));
                tindx2 = tindx+secondsToSample*sampleRate-1;
                nChans = length(OBJ.segments(1).chName);
                sc = seg_counter;
                td = 0; %total duration
                while(sc>0)
                    td = td + OBJ.segments(sc).duration * sampleRate;
                    sc = sc - 1;
                end
                if(tindx2 > td) %If the second tindx is outside of the segment, increase the segment counter.
                    seg_counter = seg_counter + 1;
                end
                data = getdataQ(OBJ, seg_counter, [tindx+secondOffset*sampleRate tindx2], desiredChannelIndices); %From event start to +5 secs from first 24 channels
                %data = data(1:secondsToSample*sampleRate); %I'm pretty sure
                %this line doesn't work as intended
                size_data = size(data);
                timeVec = max(secondOffset*sampleRate,1):secondsToSample*sampleRate-1;
                
                transposed_data = data';
                curr_signal = pop_importdata('data', transposed_data, 'srate', sampleRate, 'nbchan', length(data), 'pnts', length(transposed_data)); %Transfer to eeglab structure so we can filter with their functions
                curr_signal = pop_eegfiltnew(curr_signal, lowCutoff, highCutoff);
                curr_signal_to_plot = pop_eegfiltnew(curr_signal, bandpassLowCutoff, bandpassHighCutoff);
                % curr_signal.data = curr_signal.data(samplingRate:end);
                % hann function
                nx = max(size_data);
                na = 8;
                w = hanning(floor(nx/na));
                [pxx, f] = pwelch(curr_signal.data(1, :), w , 0 , 2.^(nextpow2(max(size_data))+1), OBJ.segments(seg_counter).samplingRate(1)); %modifying this
                
                o1pkmax = 0;
                o2pkmax = 0;
                [data, locs] = findpeaks(pxx,f); %find peaks
                [data2, locs2] = findpeaks(pxx2,f2);
                o1pks = locs(locs > lowCutoff & locs < highCutoff)'; %find peaks in the range we want
                o2pks = locs2(locs2 > lowCutoff & locs2 < highCutoff)';
                
                % Only keep the highest amplitude peak
                for i = 1:length(o1pks)
                    if (pxx(f == o1pks(i)) > pxx(f == o1pkmax))
                        o1pkmax = o1pks(i);
                    end
                end
                for i = 1:length(o2pks)
                    if (pxx(f == o2pks(i)) > pxx(f == o2pkmax))
                        o2pkmax = o2pks(i);
                    end
                end
                o1pks = o1pkmax;
                o2pks = o2pkmax;
                if (o1pks == 0)
                    manualFiles(length(manualFiles)) = [];
                    continue;
                else
                    autoEstimates = [autoEstimates rndmFName o1pks];
                end
                
                % For internal (?) use
                figure(figurecounter)
                figurecounter = figurecounter + 1;
                
                subplot(2,1,1)
                findpeaks(pxx,f); xlim([0 20]);
                % plot(f, 10*log10(pxx)); xlim([0 20]); ylim([0 30]);
                % ylabel('power spectral density   [(rad/sample)^{-1}] in
                % dB') is what is being plotted with 10*log(10)
                %hold on;
                %plot(f2, 10*log10(pxx2)); Just o1
                title(rndmFName + " Eyes Closed Event Number " + x + " " + o1pks);
                
                subplot(2,1,2)
                plot(timeVec, curr_signal.data(1,1:length(timeVec))); %Just o1
                title("Curr Signal")
                xticks([sampleRate 2*sampleRate 3*sampleRate 4*sampleRate, 5*sampleRate]) %Right now you'll have to physically change this if you increase the number of seconds
                xticklabels({'1 Second','2 Seconds','3 Seconds','4 Seconds', '5 Seconds'})
                legend(desiredChannels);
                
                figurecounter = figurecounter + 1;
                fn = "Figure " + (figurecounter-2)/2;
                set(gcf,'units','points','position',[50,100,1000,400]); %x, y, width, height
                saveas(gcf, pwd + "\" + "internal\" + fn + ".png");
                
                %For readers
                figure(figurecounter)
                
                plot(timeVec, curr_signal_to_plot.data(1,1:length(timeVec)));
                title(fn);
                xticks([sampleRate 2*sampleRate 3*sampleRate 4*sampleRate, 5*sampleRate]) %Right now you'll have to physically change this if you increase the number of seconds
                xticklabels({'1 Second','2 Seconds','3 Seconds','4 Seconds', '5 Seconds'})
                hold on;
                yl = ylim;
                if(yl(1) < -100 || yl(2) > 100)
                    ylims = [-150 150];
                elseif (yl(1) < -50 || yl(2) > 50)
                    ylims = [-100 100];
                else
                    ylims = [-50 50];
                end
                ylim(ylims)
                xlim([0 6*sampleRate])
                for i = 1:secondsToSample
                   line([i*sampleRate i*sampleRate], [ylims(1) ylims(2)], 'LineStyle', '-', 'Color', 'red');
                   hold on;
                end
                set(gcf,'units','points','position',[50,100,1000,200]); %x, y, width, height
                saveas(gcf, pwd + "\" + "readers\" + fn + ".png");
                
            catch ME
                warning(file.name+"unsuitable for random file sampling, moving on. "+ getReport(ME))
                manualFiles(length(manualFiles)) = [];
            end
        end
    end
end

%save('lastrun.mat'); %Saves workspace to file
disp("Total number of files: " + fileCounter)
disp("Total files read successfully: " + (find(filesRead=="", 1, 'first')-1))
disp("Total files read with no detected 'eyes closed' events: " + (length(noEyesClosed)-1))

%This function returns the appropriate time index given an event index
function out = getTimeIndexOfEvent(obj, eventIndex)
    starttime = obj.segments.dateOLE; % The time that the first segment of the EEG starts
    eventtime = obj.eventMarkers(eventIndex).dateOLE; % Reference event time
    diff = eventtime - starttime; %Making these variables was only necessary for my clarity
    diff = diff*3600*24; %In seconds
    arr = getSampleRate(obj, 1, 1);
    out = fix(diff*arr); %Assume all relevant channels have same sample rate... Fix to get rid of weird decimals
end

%Taken without modification from https://www.mathworks.com/matlabcentral/answers/375571-how-do-i-access-hex-values-stored-in-a-char-in-matlab
function Out = base64toHex(In)
    P  = [65:90, 97:122, 48:57, 43, 47] + 1;  % [0:9, a:z, A:Z, +, /]
    v8 = [128, 64, 32, 16, 8, 4, 2, 1];
    v6 = [32; 16; 8; 4; 2; 1];
    Table    = zeros(1, 256);
    Table(P) = 1:64;
    Data     = Table(In(:).' + 1) - 1;
    X        = rem(floor(bsxfun(@rdivide, Data, v6)), 2);
    Num      = v8 * reshape(X(1:fix(numel(X) / 8) * 8), 8, []);
    Out = sprintf('%.2x', Num);
end