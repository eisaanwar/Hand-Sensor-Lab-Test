acc% single trial
classdef trialClass
    properties
        data; %the coda files
        handDeviceData; %files from hand device if applicable (separate class for this?)
        
        massLoad; %mass the user is lifting
        participantMass;
        
        humanForcePlate;
        weightForcePlateForce;
        
        gastrocnemius; bicep; tibialisAnterior;
        
        ankle; knee; hip; shoulder; elbow; wrist; calc; MT; toe; stand; weight;
        
        timeVectorMarkers;
        timeVectorForcePlate;
        timeVectorEMG;
        
        comSystem; %x and y com
        comsSegment; %x and y com for all segments
        comWeightC; %corrected (need better name)
        weightForcePlateForceMarkerSync;
        
        segmentAngleAnkle;
        
        numberOfRepetitions;
        repetionLengths;
        maxRepetitionLength; %to set the size of the matrices
        
        weightMovementStartIndexes;
        weightMovementStartIndexesHandDevice;
        weightMimimumDisplacements;
        
        weightMinimumDistanceMarker; %should use interpolation for this.
        handDeviceDiantanceWithAngle;
        
        weightOffsetMarker;
        weightOffsetHandDevice;
    end
    
    methods
        function this = trialClass(dataLocation, handDeviceDataLocation, participantMass, massLoad, reorient) %constructor
            this.data = load(dataLocation); % load Coda data file (before cut)  'CodaFiles/LabTest10_032_2022_6_21'
            this.handDeviceData = load(handDeviceDataLocation);
            
            if(reorient) % if sideways (facing the door)
                this.data = this.reorientData(this.data); % for force plates and markers
            end
            
            this.participantMass = participantMass;
            this.massLoad = massLoad;
            
            [this.handDeviceData, this.humanForcePlate, this.weightForcePlateForce, this.gastrocnemius, this.bicep, this.tibialisAnterior, this.ankle, this.knee, this.hip, ...
                this.shoulder, this.elbow, this.wrist, this.calc, this.MT, this.toe, this.stand, this.weight] = this.cutAndFilterData(this.data, this.handDeviceData, this.massLoad);
            
            this.timeVectorMarkers = (1:length(this.ankle))/200; %to be able to plot the data in seconds (fine that this is filterned data, yes)
            this.timeVectorForcePlate = (1:length(this.humanForcePlate))/2000;
            this.timeVectorEMG = (1:length(this.gastrocnemius))/2000; % time for 1 emg sensor
            
            %calcualte CoM of segments and system %comsSegment: 1=foot, 2=lowerLeg, 3=thigh, 4=truck%Head, 5=upperArm, 6=forearm, 7=weight
            [this.comSystem, this.comsSegment, this.comWeightC, this.weightForcePlateForceMarkerSync, this.segmentAngleAnkle] = this.calcualteCoMs(this.ankle, this.calc, this.toe, this.knee, this.hip, this.shoulder, this.elbow, this.wrist, this.stand, this.weight, this.weightForcePlateForce, this.participantMass); %all 7 CoMs for x and z axis.
            
            %slice data
            [this.numberOfRepetitions, this.weightMovementStartIndexes, this.weightMovementStartIndexesHandDevice, this.weightMimimumDisplacements, this.repetionLengths, this.maxRepetitionLength] = this.sliceData(this.weight, this.handDeviceData);
            
            %hand device distance
            this.handDeviceDiantanceWithAngle = (this.handDeviceData(:,3)/1000).*cos(deg2rad(this.handDeviceData(:,4))); %/1000 for grams
            
            %Closest distance of the weight to chest, chnage to interpolate markers
            this.weightMinimumDistanceMarker = min(this.weightMimimumDisplacements);
            
            %calcaule the CoM Offset from the markers and the hand device
            this.weightOffsetMarker =  ((this.comWeightC-this.weightMinimumDistanceMarker)/1000).*(-0.1*this.weightForcePlateForceMarkerSync); % weight*distance (from the marker and force plate)
            this.weightOffsetHandDevice = (this.handDeviceData(:,2)/1000).*this.handDeviceDiantanceWithAngle;            
        end
        
        function [handDeviceData, humanForcePlate, weightForcePlateForce, gastrocnemius, bicep, tibialisAnterior, ankle, knee, hip, shoulder, elbow, wrist, calc, MT, toe, stand, weight] = ...
                cutAndFilterData(this, data, handDeviceData, massLoad)
            updateFrequencyMarker = data.Marker.RAnk.Rate; %200Hz
            updateFrequencyForcePlate = data.CoP.x1.Rate; %2000Hz
            updateFrequencyEMG = data.Analog.EMG09.Rate; %2000Hz
            
            
            %force plates centre of pressure and force
            weightForcePlateForce = this.filterForcePlate(data.Force.F2.value);
            
            %find the times the weight is lifted. (this is the data we need to use
            weightHoldingTime = find(weightForcePlateForce(:,3)<((massLoad-0.5)*-10))/updateFrequencyForcePlate; %times weight is held in seconds. Mass * 10 for weight in newtons, -1 kg for times
            weightPickupTime = weightHoldingTime(1);
            weightPutdownTime = weightHoldingTime(length(weightHoldingTime)); 
            
            markerStartIndex = round(weightPickupTime*updateFrequencyMarker); % round because force plate update frequency x10 higher
            markerEndIndex = round(weightPutdownTime*updateFrequencyMarker);
            
            emgStartIndex = int32(weightPickupTime*updateFrequencyEMG);
            engEndIndex = int32(weightPutdownTime*updateFrequencyEMG);
            
            forcePlateStartIndex = int32(weightPickupTime*updateFrequencyForcePlate);
            forcePlateEndIndex = int32(weightPutdownTime*updateFrequencyForcePlate);
            
            %find the closest find read from the hand sensor corresponding to this time,
            [handDeviceStartTimeDifference, handDeviceStartIndex] = min(abs(handDeviceData(:,1)-(weightPickupTime*1000))); % *1000 because hand hensor in ms
            [handDeviceEndTimeDifference, handDeviceEndIndex] = min(abs(handDeviceData(:,1)-(weightPutdownTime*1000)));
            
            %cut the hand data,
            handDeviceData = handDeviceData(handDeviceStartIndex:handDeviceEndIndex,:);
            %reset the time
            handDeviceData(:,1) = handDeviceData(:,1)-weightPickupTime*1000;
            
            %filterData
            humanForcePlate = this.filterForcePlate(data.CoP.x1.value(forcePlateStartIndex:forcePlateEndIndex,:));
            weightForcePlateForce = weightForcePlateForce(forcePlateStartIndex:forcePlateEndIndex,:);
            
            %EMG data
            gastrocnemius = this.envelope_emg(this.filter_emg(data.Analog.EMG09.value(:,emgStartIndex:engEndIndex)));
            bicep = this.envelope_emg(this.filter_emg(data.Analog.EMG03.value(:,emgStartIndex:engEndIndex)));
            tibialisAnterior = this.envelope_emg(this.filter_emg(data.Analog.EMG05.value(:,emgStartIndex:engEndIndex)));
            
            %markers (only get when the thing has been picked up)
            ankle = this.filterMarker(data.Marker.RAnk.value(markerStartIndex:markerEndIndex,:)); %to get sole of foot for lower extremity calulation?
            knee = this.filterMarker(data.Marker.RKnee.value(markerStartIndex:markerEndIndex,:));
            hip = this.filterMarker(data.Marker.RHip.value(markerStartIndex:markerEndIndex,:));
            shoulder = this.filterMarker(data.Marker.RShoul.value(markerStartIndex:markerEndIndex,:));
            elbow = this.filterMarker(data.Marker.RElbow.value(markerStartIndex:markerEndIndex,:));
            wrist = this.filterMarker(data.Marker.RWrist.value(markerStartIndex:markerEndIndex,:));
            calc = this.filterMarker(data.Marker.RCalc.value(markerStartIndex:markerEndIndex,:));
            MT = this.filterMarker(data.Marker.R5th.value(markerStartIndex:markerEndIndex,:));
            toe = this.filterMarker(data.Marker.RToe.value(markerStartIndex:markerEndIndex,:));
            stand = this.filterMarker(data.Marker.Stand.value(markerStartIndex:markerEndIndex,:));
            weight = this.filterMarker(data.Marker.Weight.value(markerStartIndex:markerEndIndex,:));
        end
        
        function [numberOfRepetitions, weightMovementStartIndexes, weightMovementStartIndexesHandDevice, weightMimimumDisplacements, repetionLengths, maxRepetitionLength] = ...
                sliceData(this, weight, handDeviceData)
            % can also make sure that the points are last valey to next peak (works
            % with muliple peaks in the same place
            [weightMimimumDisplacements, weightMovementStartIndexes] = findpeaks(-weight(:,1));
            weightMimimumDisplacements = -weightMimimumDisplacements; % negative becuase inverted then found peaks
            [weightMaximumDisplacements, weightMovementPeakIndexes] = findpeaks(weight(:,1)); %peaks are not being used yet, only to see number of repetiotns
            
            %calcualte what max and min represent the movements
            if weightMovementPeakIndexes(1) < weightMovementStartIndexes(1) % there shouldnt be an end point with no start at the beginning
                weightMovementPeakIndexes(1) = [];
            end
            
            if weightMovementPeakIndexes(size(weightMovementPeakIndexes)) > weightMovementStartIndexes(size(weightMovementStartIndexes)) % last point should be min (if a peak then moved back as put down)
                weightMovementPeakIndexes(size(weightMovementPeakIndexes)) = [];
            end
            
            if weightMovementStartIndexes(end) > weightMovementPeakIndexes(end) %calcualte the number of repetitions
                numberOfRepetitions = length(weightMovementStartIndexes)-1;
            else
                numberOfRepetitions = length(weightMovementStartIndexes);
            end
            
            weightMovementStartTimes = weightMovementStartIndexes/200;
            weightMovementPeakTimes = weightMovementPeakIndexes/200; %this is the peak point, could be the end point
            
            % hand sensor start indexes
            weightMovementStartIndexesHandDevice = NaN(numberOfRepetitions,1);
            for i = 1:numberOfRepetitions+1 %add one becaues need the last number with is not a full repetition
                [handDeviceTimeDifference, weightMovementStartIndexeHandDevice] = min(abs(handDeviceData(:,1)-(weightMovementStartTimes(i)*1000))); % *1000 because hand hensor in ms, find closest point in data
                weightMovementStartIndexesHandDevice(i) = weightMovementStartIndexeHandDevice;
            end
            
            % find the longest repetion to use as the matrix size
            maxRepetitionLength = 0;
            for i = 1:numberOfRepetitions
                repetionLength = weightMovementStartIndexes(i+1)-weightMovementStartIndexes(i);
                if maxRepetitionLength < repetionLength
                    maxRepetitionLength = repetionLength;
                end
            end
            
            % find the longest repetion to use as the matrix sizes for the hand sensor
            maxRepetitionLengthHandSensor = 0;
            for i = 1:numberOfRepetitions
                repetionLength = weightMovementStartIndexesHandDevice(i+1)-weightMovementStartIndexesHandDevice(i) ;
                if maxRepetitionLength < repetionLength
                    maxRepetitionLengthHandSensor = repetionLength;
                end
            end
            
            repetionLengths = NaN(numberOfRepetitions,1);
            for i = 1:numberOfRepetitions
                repetionLengths(i) = weightMovementStartIndexes(i+1)-weightMovementStartIndexes(i);
            end
            
            repetionLengthsHandSensor = NaN(numberOfRepetitions,1);
            for i = 1:numberOfRepetitions
                repetionLengthsHandSensor(i) = weightMovementStartIndexesHandDevice(i+1)-weightMovementStartIndexesHandDevice(i);
            end
        end
        
        function data = reorientData(this, data)
            
            %force plates
            %swap axis
            temp = data.CoP.x1.value(:,1);
            data.CoP.x1.value(:,1) = -data.CoP.x1.value(:,2);
            data.CoP.x1.value(:,2) = temp;
            
            %ankle
            temp = data.Marker.RAnk.value(:,1);
            data.Marker.RAnk.value(:,1) = -data.Marker.RAnk.value(:,2); % swap the x and y axis
            data.Marker.RAnk.value(:,2) = temp;
            
            %knee
            temp = data.Marker.RKnee.value(:,1);
            data.Marker.RKnee.value(:,1) = -data.Marker.RKnee.value(:,2); % swap the x and y axis
            data.Marker.RKnee.value(:,2) = temp;
            
            %hip
            temp = data.Marker.RHip.value(:,1);
            data.Marker.RHip.value(:,1) = -data.Marker.RHip.value(:,2); % swap the x and y axis
            data.Marker.RHip.value(:,2) = temp;
            
            %Shoulder
            temp = data.Marker.RHip.value(:,1);
            data.Marker.RShoul.value(:,1) = -data.Marker.RShoul.value(:,2); % swap the x and y axis
            data.Marker.RShoul.value(:,2) = temp;
            
            %elbow
            temp = data.Marker.RElbow.value(:,1);
            data.Marker.RElbow.value(:,1) = -data.Marker.RElbow.value(:,2); % swap the x and y axis
            data.Marker.RElbow.value(:,2) = temp;
            
            %wrist
            temp = data.Marker.RWrist.value(:,1);
            data.Marker.RWrist.value(:,1) = -data.Marker.RWrist.value(:,2); % swap the x and y axis
            data.Marker.RWrist.value(:,2) = temp;
            
            %calc
            temp = data.Marker.RCalc.value(:,1);
            data.Marker.RCalc.value(:,1) = -data.Marker.RCalc.value(:,2); % swap the x and y axis
            data.Marker.RCalc.value(:,2) = temp;
            
            %MT
            temp = data.Marker.R5th.value(:,1);
            data.Marker.R5th.value(:,1) = -data.Marker.R5th.value(:,2); % swap the x and y axis
            data.Marker.R5th.value(:,2) = temp;
            
            %RToe
            temp = data.Marker.RToe.value(:,1);
            data.Marker.RToe.value(:,1) = -data.Marker.RToe.value(:,2); % swap the x and y axis
            data.Marker.RToe.value(:,2) = temp;
            
            %Stand
            temp = data.Marker.Stand.value(:,1);
            data.Marker.Stand.value(:,1) = -data.Marker.Stand.value(:,2); % swap the x and y axis
            data.Marker.Stand.value(:,2) = temp;
            
            %Weight
            temp = data.Marker.Weight.value(:,1);
            data.Marker.Weight.value(:,1) = -data.Marker.Weight.value(:,2); % swap the x and y axis
            data.Marker.Weight.value(:,2) = temp;
            
            
        end        

        function [comSystem, comsSegment, comWeightC, weightForcePlateForceMarkerSync, segmentAngleAnkle] = ...
                calcualteCoMs(this, ankle, calc, toe, knee, hip, shoulder, elbow, wrist, stand, weight, weightForcePlateForce, participantMass)
            %calcualte CoM of segments and system
            %1=foot, 2=lowerLeg, 3=thigh, 4=truck%Head, 5=upperArm, 6=forearm&Hand,
            %7=weight, 8=stand
            comsSegment = NaN(length(ankle(:,1)),2,8); %all 7 CoMs for x and z axis.
            
            %foot CoM
            footLength = -(calc(:,1) - toe(:,1)); %what if the foot is tilted
            footHeight = ankle(:,3); %what if the ankle is lifted (use the other angles)
            comsSegment(:,1,1) = calc(:,1) + (0.4485*footLength); % the x postion is just the ankle position and the y postion is between the ankle and floor (not exactly)
            comsSegment(:,2,1) = ankle(:,3) - (0.4622*footHeight); %vertical component
            
            %lower leg CoM
            segmentAngleAnkle = atan2((knee(:,1)-ankle(:,1)),(knee(:,3)-ankle(:,3))); % ankle angles
            lowerlegLength = sqrt((knee(:,1)-ankle(:,1)).^2 + (knee(:,3)-ankle(:,3)).^2); %is it inefficient to keep calling the variable?
            comsSegment(:,1,2) = ankle(:,1) + (0.6295*lowerlegLength).*sin(segmentAngleAnkle); %xCoM
            comsSegment(:,2,2) = ankle(:,3) + (0.6295*lowerlegLength).*cos(segmentAngleAnkle); %yCoM
            
            %thigh CoM
            kneeAngle = atan2((hip(:,1)-knee(:,1)),(hip(:,3)-knee(:,3))); % knee angles
            thighLength = sqrt((hip(:,1)-knee(:,1)).^2 + (hip(:,3)-knee(:,3)).^2);  % thigh length
            comsSegment(:,1,3) = knee(:,1) + (0.6281*thighLength).*sin(kneeAngle); % thigh xCoM
            comsSegment(:,2,3) = knee(:,3) + (0.6281*thighLength).*cos(kneeAngle); % thigh yCoM
            
            %truckAndHead CoM
            hipAngle = atan2((shoulder(:,1)-hip(:,1)),(shoulder(:,3)-hip(:,3))); % knee angles
            truckAndHeadLength = sqrt((shoulder(:,1)-hip(:,1)).^2 + (shoulder(:,3)-hip(:,3)).^2);  % truck and head length, marker goes to shoulder but the promial and distal is different
            comsSegment(:,1,4) = hip(:,1) + (0.5021*truckAndHeadLength).*sin(hipAngle); % truckAndHead xCoM, used different distal
            comsSegment(:,2,4) = hip(:,3) + (0.5021*truckAndHeadLength).*cos(hipAngle); % truckAndHead yCoM
            
            %upperArm CoM
            upperArmAngle = atan2((shoulder(:,1)-elbow(:,1)),(shoulder(:,3)-elbow(:,3))); % upper arm angles
            upperArmLength = sqrt((shoulder(:,1)-elbow(:,1)).^2 + (shoulder(:,3)-elbow(:,3)).^2);  % upper arm length
            comsSegment(:,1,5) = shoulder(:,1) + (0.5130*upperArmLength).*-sin(upperArmAngle); % xCoM
            comsSegment(:,2,5) = shoulder(:,3) + (0.5130*upperArmLength).*-cos(upperArmAngle); % yCoM
            
            %forearm CoM
            forarmAngle = atan2((elbow(:,1)-wrist(:,1)),(elbow(:,3)-wrist(:,3))); % forearm angles
            forearmLength = sqrt((elbow(:,1)-wrist(:,1)).^2 + (elbow(:,3)-wrist(:,3)).^2);  % forearm length
            comsSegment(:,1,6) = elbow(:,1) + (0.3896*forearmLength).*-sin(forarmAngle); % xCoM
            comsSegment(:,2,6) = elbow(:,3) + (0.3896*forearmLength).*-cos(forarmAngle); % yCoM
            
            %weight CoM (need system to calculate the COM of weight (need 2 markers)
            comsSegment(:,1,7) = weight(:,1); %need for for this (use stand but not visable)
            comsSegment(:,2,7) = weight(:,3);
            
            %1=foot, 2=lowerLeg, 3=thigh, 4=truck%Head, 5=upperArm, 6=forearm&Hand,
            massProportions = [2*0.0147; 2*0.0435; 2*0.1027; 0.5801; 2*0.0263; 2*0.0227;]*participantMass;
            
            %CHANGE NAMESSSSS!!! (to sync markers with the force plate)
            counts = floor(length(weightForcePlateForce(:,3))/10);
            comWeightC = comsSegment(1:counts,1,7);
            weightForcePlateForceMarkerSync = weightForcePlateForce(1:10:counts*10,3); %use 1 in every 10 readings
            
            %calculate the CoM
            comSystem(:,1) = comsSegment(1:counts,1,1).*massProportions(1) + comsSegment(1:counts,1,2).*massProportions(2) ...
                + comsSegment(1:counts,1,3).*massProportions(3) + comsSegment(1:counts,1,4).*massProportions(4) ...
                + comsSegment(1:counts,1,5).*massProportions(5) + comsSegment(1:counts,1,6).*massProportions(6) ...
                + comWeightC.*(-0.1*weightForcePlateForceMarkerSync); %(need to make all shorter then)
            
            comSystem(:,1) = comSystem(:,1)./(participantMass+(0.1*-weightForcePlateForceMarkerSync));
            % divid by total mass (maybe add a bit for the hand sensor (stand is)
            % included, but not the backplate
        end
        
        function plotComs(this, calc, ankle, MT, toe, knee, hip, shoulder, elbow, wrist, comsSegment)
            %lower body CoMs
            figure('name','CoM Lower Body Segements');
            tiledlayout(3,2)
            
            %foot com
            nexttile
            plot(comsSegment(:,1,1));
            hold on;
            plot(calc(:,1));
            plot(ankle(:,1));
            plot(MT(:,1));
            plot(toe(:,1));
            legend('xCoM Foot','Calc X', 'Ankle X', 'MT X', 'Toe X')
            title('xCoM Foot')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,1));
            hold on;
            plot(calc(:,3));
            plot(ankle(:,3));
            plot(MT(:,3));
            plot(toe(:,3));
            legend('yCoM Foot','Calc Y', 'Ankle Y', 'MT Y', 'Toe X')
            title('yCoM Foot')
            xlabel('Counts')
            ylabel('Distance (mm)')
            
            
            %lower leg CoM
            nexttile
            plot(comsSegment(:,1,2));
            hold on;
            plot(ankle(:,1));
            plot(knee(:,1));
            legend('xCoM Lower Leg','Ankle X', 'Knee X')
            title('xCoM Lower Leg')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,2));
            hold on;
            plot(ankle(:,3));
            plot(knee(:,3));
            legend('yCoM Lower Leg','Ankle Y', 'Knee Y')
            title('yCoM Lower Leg')
            xlabel('Counts')
            ylabel('Distance (mm)')
            
            %thigh CoM
            nexttile
            plot(comsSegment(:,1,3));
            hold on;
            plot(knee(:,1));
            plot(hip(:,1));
            legend('xCoM Thigh','Knee X', 'Hip X')
            title('xCoM Thigh')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,3));
            hold on;
            plot(knee(:,3));
            plot(hip(:,3));
            legend('yCoM Thigh','Knee Y', 'Hip Y')
            title('yCoM Thigh')
            xlabel('Counts')
            ylabel('Distance (mm)')
            
            %Upper body CoMs
            figure('name','Upper Body Segments');
            tiledlayout(3,2)
            
            %trunkAndHead CoM
            nexttile
            plot(comsSegment(:,1,4));
            hold on;
            plot(hip(:,1));
            plot(shoulder(:,1));
            legend('xCoM TrunkAndHead','Hip X', 'Shoulder X')
            title('xCoM TrunkAndHead')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,4));
            hold on;
            plot(hip(:,3));
            plot(shoulder(:,3));
            legend('yCoM TrunkAndHead','Hip Y', 'Shoulder Y')
            title('yCoM TrunkAndHead')
            xlabel('Counts')
            ylabel('Distance (mm)')
            
            %upperArm CoM
            nexttile
            plot(comsSegment(:,1,5));
            hold on;
            plot(shoulder(:,1));
            plot(elbow(:,1));
            legend('xCoM UpperArm','Shoulder X', 'Elbow X')
            title('xCoM UpperArm')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,5));
            hold on;
            plot(shoulder(:,3));
            plot(elbow(:,3));
            legend('yCoM TrunkAndHead','Shoulder Y', 'Elbow Y')
            title('yCoM TrunkAndHead')
            xlabel('Counts')
            ylabel('Distance (mm)')
            
            %forearm CoM
            nexttile
            plot(comsSegment(:,1,6));
            hold on;
            plot(elbow(:,1));
            plot(wrist(:,1));
            legend('xCoM Forearm','Elbow X', 'Wrist X')
            title('xCoM Forearm')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,6));
            hold on;
            plot(elbow(:,3));
            plot(wrist(:,3));
            legend('yCoM Forearm','Elbow Y', 'Wrist Y')
            title('yCoM Forearm')
            xlabel('Counts')
            ylabel('Distance (mm)')
            
            %weight CoMs
            figure('name','CoM Weights');
            tiledlayout(2,1)
            nexttile
            plot(comsSegment(:,1,7));
            hold on;
            legend('xCoM Weight')
            title('xCoM Weight')
            xlabel('Counts')
            ylabel('Distance (mm)')
            nexttile
            plot(comsSegment(:,2,7));
            hold on;
            legend('yCoM Weight')
            title('yCoM Weight')
            xlabel('Counts')
            ylabel('Distance (mm)')
        end
        
        function plotWeightDistance(this, comWeight, timeVectorMarkers, weightMinimumDistanceMarker, handDeviceData, handDeviceDiantanceWithAngle)
            figure('name','Weight Distance');
            tiledlayout(3,1)
            
            nexttile
            plot(timeVectorMarkers,((comWeight(:,1)-weightMinimumDistanceMarker)/1000)); %/1000 for meters
            ylim([0 0.75])
            title('Marker Distance')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile;
            plot(handDeviceData(:,1)/1000,handDeviceData(:,3));
            ylim([0 0.75])
            title('Hand Device Without Angle')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile;
            plot(handDeviceData(:,1)/1000,handDeviceDiantanceWithAngle);
            ylim([0 0.75])
            title('Hand Device With Angle')
            xlabel('Time (s)')
            ylabel('Distance (m)')
        end
        
        function plotComOffset(this, numberOfRepetitions, weightOffsetHandDevice, weightMovementStartIndexesHandDevice, ...
                handDeviceData, timeVectorMarkers, weightOffsetMarker)
            
            figure('name','CoM Offset');
            tiledlayout(2,1)
            nexttile
            plot(timeVectorMarkers(1:length(weightOffsetMarker)), weightOffsetMarker);
            ylim([0 2.5])
            title('CoM Offset Marker')
            xlabel('Time (s)')
            ylabel('Offset (Nm)')
            
            nexttile;
            plot(handDeviceData(:,1)/1000, weightOffsetHandDevice);
            ylim([0 2.5])
            title('CoM Offset Hand Device')
            xlabel('Time (s)')
            ylabel('Offset (Nm)')
            
            
            figure('name','Sliced CoM Offset');
            tiledlayout(2,1)
            nexttile
            for i = 1:numberOfRepetitions
                plot(handDeviceData(1:weightMovementStartIndexesHandDevice(i+1)-weightMovementStartIndexesHandDevice(i)+1,1),... % time
                    weightOffsetHandDevice(weightMovementStartIndexesHandDevice(i):weightMovementStartIndexesHandDevice(i+1))) % com distance
                hold on;
            end
            title('Hand Sensor CoM Horizontal Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile
            for i = 1:numberOfRepetitions
                plot(handDeviceData(1:weightMovementStartIndexesHandDevice(i+1)-weightMovementStartIndexesHandDevice(i)+1,1),... % time
                    weightOffsetHandDevice(weightMovementStartIndexesHandDevice(i):weightMovementStartIndexesHandDevice(i+1))) % com distance
                hold on;
            end
            title('Marker CoM Horizontal Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
        end
        
        function plotSlicedMovementData(this, timeVectorMarkers, timeVectorForcePlate, maxRepetitionLength, numberOfRepetitions, ...
                repetionLengths, humanForcePlate, xComTotalSystem, weight, weightMovementStartIndexes, segmentAngleAnkle)
            
            figure('name','Sliced Marker Data');
            tiledlayout(2,2)
            
            nexttile
            weightHorizontalMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                weightHorizontalMatrix(1:repetionLengths(i)+1,i) = weight(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1),1);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), weightHorizontalMatrix);
            title('Weight Horizontal Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile
            %segmentAngleAnkle = transpose(segmentAngleAnkle);
            segmentAngleAnkle = rad2deg(segmentAngleAnkle);
            ankleAngleMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                ankleAngleMatrix(1:repetionLengths(i)+1,i) = segmentAngleAnkle(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1),1);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), ankleAngleMatrix);
            title('Ankle Angle for Each Repetition')
            xlabel('Time (s)')
            ylabel('Angle (degrees)')
            
            nexttile
            %xComTotalSystem = transpose(xComTotalSystem);
            systemComMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                systemComMatrix(1:repetionLengths(i)+1,i) = xComTotalSystem(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1));
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), systemComMatrix);
            title('Horizontal CoM for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile
            %need index multipler
            %xComTotalSystem = transpose(xComTotalSystem);
            humanForcePlateMatrix = NaN(maxRepetitionLength*10,numberOfRepetitions); %force plate refreses 10x amount
            for i = 1:numberOfRepetitions
                humanForcePlateMatrix(1:repetionLengths(i)*10+1,i) = humanForcePlate(weightMovementStartIndexes(i)*10:weightMovementStartIndexes(i+1)*10,1);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorForcePlate(1:maxRepetitionLength*10+1), humanForcePlateMatrix);
            title('Horizontal CoP for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
        end
        
        function plotSlicedEmg(this, timeVectorMarkers, timeVectorForcePlate, repetionLengths, maxRepetitionLength, numberOfRepetitions, ...
                weightMovementStartIndexes, weight, gastrocnemius, tibialisAnterior, bicep, segmentAngleAnkle)
            
            figure('name','Sliced EMG Data');
            tiledlayout(5,1)
            
            nexttile
            %segmentAngleAnkle = transpose(segmentAngleAnkle);
            segmentAngleAnkle = rad2deg(segmentAngleAnkle);
            ankleAngleMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                ankleAngleMatrix(1:repetionLengths(i)+1,i) = segmentAngleAnkle(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1),1);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), ankleAngleMatrix);
            title('Ankle Angle for Each Repetition')
            xlabel('Time (s)')
            ylabel('Angle (degrees)')
            
            nexttile
            %need index multipler
            %xComTotalSystem = transpose(xComTotalSystem);
            gastrocnemiusMatrix = NaN(maxRepetitionLength*10,numberOfRepetitions); %force plate refreses 10x amount
            for i = 1:numberOfRepetitions
                gastrocnemiusMatrix(1:repetionLengths(i)*10+1,i) = gastrocnemius(weightMovementStartIndexes(i)*10:weightMovementStartIndexes(i+1)*10);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorForcePlate(1:maxRepetitionLength*10+1), gastrocnemiusMatrix);
            title('Gastrocnemius Actvity for Each Repetition')
            xlabel('Time (s)')
            ylabel('Voltage (mV)')
            
            nexttile
            %need index multipler
            %xComTotalSystem = transpose(xComTotalSystem);
            tibialisAnteriorMatrix = NaN(maxRepetitionLength*10,numberOfRepetitions); %force plate refreses 10x amount
            for i = 1:numberOfRepetitions
                tibialisAnteriorMatrix(1:repetionLengths(i)*10+1,i) = tibialisAnterior(weightMovementStartIndexes(i)*10:weightMovementStartIndexes(i+1)*10);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorForcePlate(1:maxRepetitionLength*10+1), tibialisAnteriorMatrix);
            title('Tibialis Anterior Actvity for Each Repetition')
            xlabel('Time (s)')
            ylabel('Voltage (mV)')
            
            nexttile
            weightMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                weightMatrix(1:repetionLengths(i)+1,i) = weight(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1),1);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), weightMatrix);
            title('Weight Horizontal Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile
            %need index multipler
            %xComTotalSystem = transpose(xComTotalSystem);
            bicepMatrix = NaN(maxRepetitionLength*10,numberOfRepetitions); %force plate refreses 10x amount
            for i = 1:numberOfRepetitions
                bicepMatrix(1:repetionLengths(i)*10+1,i) = bicep(weightMovementStartIndexes(i)*10:weightMovementStartIndexes(i+1)*10);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorForcePlate(1:maxRepetitionLength*10+1), bicepMatrix);
            title('Bicep Actvity for Each Repetition')
            xlabel('Time (s)')
            ylabel('Voltage (mV)')
        end
        
        function plotSummaryData(this, calc, MT, toe, timeVectorMarkers, timeVectorForcePlate, comSystemX, weight, weightMinimumDistanceMarker, segmentAngleAnkle,...
                humanForcePlate, weightForcePlateForce, handDeviceData, handDeviceDiantanceWithAngle, weightOffsetMarker, weightOffsetHandDevice)
            
            figure('name','Summary Data');
            tiledlayout(3,2) % 2 by 3 layout
            
            %ankle angle
            nexttile
            segmentAngleAnkle = rad2deg(segmentAngleAnkle);
            %y = segmentAngleAnkle;
            %x = timeAxis;
            plot(timeVectorMarkers, segmentAngleAnkle)
            title('Ankle Angle')
            xlabel('Time (s)')
            ylabel('Angle (degrees)')
            
            %overall CoM
            nexttile
            plot(timeVectorMarkers(1:length(comSystemX)), comSystemX, 'b')
            hold on;
            plot(timeVectorForcePlate, humanForcePlate(:,1), 'r')
            plot(timeVectorMarkers(1:length(calc(:,1))), calc(:,1), 'g') %back of base of support (foot)
            plot(timeVectorMarkers(1:length(toe(:,1))), toe(:,1), 'm') %front of base of support (foot)
            plot(timeVectorMarkers(1:length(MT(:,1))), MT(:,1), 'c')
            %plot(timeVectorMarkers, xComTotalSystem, 'r')
            %plot(timeAxisMarkers, xComBody, 'b')
            title('Horizontal Centre of Mass (CoM)')
            xlabel('Time (s)')
            ylabel('Distance (mm)')
            legend('CoM (Markers)','CoM (Force Plate)','Calc','Toe')
            
            %force plate
            nexttile
            plot(timeVectorForcePlate, humanForcePlate(:,1), 'r')
            title('Force Plate CoP')
            xlabel('Time (s)')
            ylabel('Distance (mm)')
            hold on
            plot(timeVectorForcePlate, humanForcePlate(:,2), 'b')
            legend('Anterior/Posterior CoP','Lateral/Medial CoP')
            
            nexttile
            plot(handDeviceData(:,1)/1000, handDeviceDiantanceWithAngle, 'r')
            title('Weight Horizontal Displacement')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            hold on;
            plot(timeVectorMarkers, (weight(:,1))/1000, 'b') % -weightMinimumDistanceMarker
            legend('Hand Sensor','Marker')
            
            nexttile
            %plot(timeAxisHandDevice, handDeviceLoad)
            plot(handDeviceData(:,1)/1000,handDeviceData(:,2))
            title('Hand Sensor Load')
            xlabel('Time ')
            ylabel('Load (g)')
            hold on;
            plot(timeVectorForcePlate, weightForcePlateForce(:,3)*-100, 'g')
            legend('Hand Sensor','Weight Force Plate')
            
            nexttile
            plot(timeVectorMarkers(1:length(weightOffsetMarker)), weightOffsetMarker);
            hold on;
            plot(handDeviceData(:,1)/1000, weightOffsetHandDevice);
            ylim([0 1.5])
            title('CoM Offset Marker')
            xlabel('Time (s)')
            ylabel('Offset (Nm)')
            legend('Marker CoM Offset','Hand Sensor CoM Offset')
        end
        
        function plotMarkerFpCom(this, calc, MT, toe, timeVectorMarkers, timeVectorForcePlate, comSystemX, humanForcePlate)
            plot(timeVectorMarkers(1:length(comSystemX)), comSystemX, 'b')
            hold on;
            plot(timeVectorForcePlate, humanForcePlate(:,1), 'r')
            plot(timeVectorMarkers(1:length(calc(:,1))), calc(:,1), 'g') %back of base of support (foot)
            plot(timeVectorMarkers(1:length(toe(:,1))), toe(:,1), 'm') %front of base of support (foot)
            plot(timeVectorMarkers(1:length(MT(:,1))), MT(:,1), 'c')
            %plot(timeVectorMarkers, xComTotalSystem, 'r')
            %plot(timeAxisMarkers, xComBody, 'b')
            %title('Horizontal Centre of Mass (CoM)')
            xlabel('Time (s)')
            ylabel('Distance (mm)')
            legend('CoM (Markers)','CoM (Force Plate)','Calc','Toe')
        end
        
        function plotEmgOverall(this, gastrocnemius, tibialisAnterior, bicep, timeVectorEMG, segmentAngleAnkle, timeVectorMarkers, ...
                handDeviceDataDistance, handDeviceDataTime)
            
            figure('name','EMG Overall');
            subplot(5,1,1)
            plot(timeVectorMarkers, segmentAngleAnkle, 'r')
            title('Ankle Angle')
            xlabel('Time (s)')
            ylabel('Angle (degrees)')
            
            subplot(5,1,2)
            plot(timeVectorEMG, gastrocnemius , 'b')
            xlabel('Time (s)')
            ylabel('Voltage (mV)')
            title('Gastrocnemius')
            
            subplot(5,1,3)
            plot(timeVectorEMG,tibialisAnterior, 'b')
            xlabel('Time (s)')
            ylabel('Voltage (mV)')
            title('Tibialis Anterior')
            
            subplot(5,1,4)
            plot(handDeviceDataTime/1000,handDeviceDataDistance)
            title('Hand Sensor Distance')
            xlabel('Time ')
            ylabel('Distance (mm)')
            
            subplot(5,1,5)
            plot(timeVectorEMG,bicep, 'b')
            xlabel('Time (s)')
            ylabel('Voltage (mV)')
            title('Bicep')
        end
        
        function plotSlicedHandSensorData(this, maxRepetitionLength, numberOfRepetitions, handDeviceData, weightMovementStartIndexes, ...
                weightMovementStartIndexesHandDevice, timeVectorMarkers, weight, repetionLengths)
            
            figure('name','Sliced Hand Sensor Data');
            tiledlayout(2,2)
            
            nexttile
            %handSensorDistanceMatrix = NaN(maxRepetitionLengthHandSensor,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                %handSensorDistanceMatrix(1:repetionLengthsHandSensor(i)+1,i) = weight(weightMovementStartIndexesHandDevice(i):weightMovementStartIndexesHandDevice(i+1),1);
                plot(handDeviceData(weightMovementStartIndexesHandDevice(i)-weightMovementStartIndexesHandDevice(i)+1:weightMovementStartIndexesHandDevice(i+1)-weightMovementStartIndexesHandDevice(i)+1,1),... %+1 otherwise looks for 0
                    handDeviceData(weightMovementStartIndexesHandDevice(i):weightMovementStartIndexesHandDevice(i+1),3))
                hold on;
            end %why plus 1? (becuase starting from 1)
            %plot(handDeviceData(:,1)/1000(1:maxRepetitionLengthHandSensor+1), handSensorDistanceMatrix);
            title('Hand Sensor Horizontal Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile
            weightMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                weightMatrix(1:repetionLengths(i)+1,i) = weight(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1),1);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), weightMatrix);
            title('Weight Marker Horizontal Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
            
            nexttile
            %handSensorWeightMatrix = NaN(maxRepetitionLengthHandSensor,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                %handSensorDistanceMatrix(1:repetionLengthsHandSensor(i)+1,i) = weight(weightMovementStartIndexesHandDevice(i):weightMovementStartIndexesHandDevice(i+1),1);
                plot(handDeviceData(weightMovementStartIndexesHandDevice(i)-weightMovementStartIndexesHandDevice(i)+1:weightMovementStartIndexesHandDevice(i+1)-weightMovementStartIndexesHandDevice(i)+1,1),... %+1 otherwise looks for 0
                    handDeviceData(weightMovementStartIndexesHandDevice(i):weightMovementStartIndexesHandDevice(i+1),2))
                hold on;
            end %why plus 1? (becuase starting from 1)
            %plot(handDeviceData(:,1)/1000(1:maxRepetitionLengthHandSensor+1), handSensorForceMatrix);
            title('Hand Sensor Force for Each Repetition')
            xlabel('Time (s)')
            ylabel('Weight (kg)')
            
            nexttile
            weightVerticalMatrix = NaN(maxRepetitionLength,numberOfRepetitions);
            for i = 1:numberOfRepetitions
                weightVerticalMatrix(1:repetionLengths(i)+1,i) = weight(weightMovementStartIndexes(i):weightMovementStartIndexes(i+1),3);
            end %why plus 1? (becuase starting from 1)
            plot(timeVectorMarkers(1:maxRepetitionLength+1), weightVerticalMatrix);
            title('Weight Marker Vertical Displacement for Each Repetition')
            xlabel('Time (s)')
            ylabel('Distance (m)')
        end
        
        function plotWeightmovements(this, weight, timeVectorMarkers)
            figure('name','Weight X Movements');
            tiledlayout(4,1)
            
            nexttile
            %plot(timeVectorMarkers,(weight(:,1)));
            findpeaks(weight(:,1));
            title('Hand Sensor Distance')
            xlabel('Time ')
            ylabel('Distance (mm)')
            
            nexttile
            findpeaks(-weight(:,1));
            title('Inverted Hand Sensor Distance')
            xlabel('Time ')
            ylabel('Distance (-mm)')
            
            nexttile
            weightMarkerXgradient = gradient(weight(:,1));
            plot(timeVectorMarkers,weightMarkerXgradient);
            title('Hand Sensor Speed')
            xlabel('Time ')
            ylabel('Speed (mm/s)')
            
            nexttile
            plot(timeVectorMarkers,gradient(weightMarkerXgradient));
            title('Hand Sensor Acceleration')
            xlabel('Time ')
            ylabel('Acceleratin (mm/s^2)')
            
        end
        
        function emg = filter_emg(this, emg_raw)
            %global f_s;
            f_s = 2000;
            %global time;
            % Check frequency transform for peaks (noise)
            % 2-4 order butterworth
            % highpass -> rectify -> lowpass
            % 50 or 100 -> abs -> 3-4 Hz
            
            % 1. De-mean the raw signal (subtract out the mean)
            emg_demean = emg_raw - mean(emg_raw);
            % 2. High pass filter: [50/(f_s/2) -> 250 Hz = 1/(50/(f_s/2))]
            %[b, a] = butter(4, 50/(f_s/2), 'high');
            [b, a] = butter(4, 0.1, 'high'); %Eisa %200Hz high pass (20 to 500, linear envolope
            emg_highpass = filtfilt(b , a, emg_demean');
            % 3. Apply a Notch Filter at 60 Hz
            % NONE NECESSARY HERE (no spike in FFT)
            emg_notch = emg_highpass;
            
            emg = emg_notch';
            
            %   if show_plots
            %       emg_raw_f = fftshift(fft(emg_demean));
            %       dF = f_s/length(time);
            %       f_t = -f_s/2:dF:f_s/2-dF;
            %       figure; subplot(2,1,1);
            %       plot(f_t, abs(emg_raw_f)/length(time));
            %       title('Raw EMG frequency spectrum');
            %       xlim([0, f_s/2])
            %       subplot(2,1,2);
            %       emg_f = fftshift(fft(emg));
            %       plot(f_t, abs(emg_f)/length(time));
            %       xlim([0.01, 5])
            %       xlabel('Frequency (Hz)');
            %   end
        end
        
        function emg = envelope_emg(this, emg_filt)
            %global f_s;
            %global time;
            f_s = 2000;
            %time = emgAxisMarkers
            % Rectify and envelope EMG
            % 4. Rectify (take absolute value)
            emg_rect = abs(emg_filt);
            % 5. Apply an FIR LowPass filter with cut-off freq = 10 Hz (To get Linear Envelope)
            %[2/(f_s/2) -> 10 Hz = 1/(2/(f_s/2))]
            [b, a] = butter(4, 2/(f_s/2));
            emg = filtfilt(b, a, emg_rect);
        end
        
        function marker = filterMarker(this, markerData)
            [b,a] = butter(2,0.01); %0.01 is 10Hz
            marker = filtfilt(b,a,markerData);
        end
        
        function forcePlate = filterForcePlate(this, forcePlateData)
            [b,a] = butter(4,0.01);
            forcePlate = filtfilt(b,a,forcePlateData);
        end
        
        function drawAnimation(this, ankle, knee, shoulder, elbow, wrist, calc, MT, stand, weight, xComFoot, xCoMFoot)
            for i = 1:length(ankle)
                clf
                if mod(i,2) == 1
                    continue;
                end
                
                %plot(timeAxisMarkers(i), xCoMh(:,i),);
                
                %plot marker points
                
                plot(ankle(i,1), ankle(i,3), '.');
                text(ankle(i,1),ankle(i,3),'Ankle','VerticalAlignment','bottom','HorizontalAlignment','right')
                axis([0 2000 0 2000]);
                set(gcf,'position',[100,100,750,700])
                
                hold on
                plot(knee(i,1), knee(i,3), '.');
                text(knee(i,1),knee(i,3),'Knee','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(hip(i,1), hip(i,3), '.');
                text(hip(i,1),hip(i,3),'Hip','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(shoulder(i,1), shoulder(i,3), '.');
                text(shoulder(i,1),shoulder(i,3),'Shoulder','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(elbow(i,1), elbow(i,3), '.');
                text(elbow(i,1),elbow(i,3),'Elbow','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(wrist(i,1), wrist(i,3), '.');
                text(wrist(i,1),wrist(i,3),'Wrist','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(calc(i,1), calc(i,3), '.');
                text(calc(i,1),calc(i,3),'Calc','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(MT(i,1), MT(i,3), '.');
                text(MT(i,1),MT(i,3),'MT','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                plot(stand(i,1), stand(i,3), '.');
                text(stand(i,1),stand(i,3),'Stand','VerticalAlignment','bottom','HorizontalAlignment','left')
                
                plot(weight(i,1), weight(i,3), '.');
                text(weight(i,1),weight(i,3),'Weight','VerticalAlignment','bottom','HorizontalAlignment','left')
                
                %plot CoM points
                plot(xComFoot(i,1), xCoMFoot(i,3), '.');
                text(knee(i,1),knee(i,3),'xCoMFoot','VerticalAlignment','bottom','HorizontalAlignment','right')
                
                
                drawnow
                
                i = i + 100; %draw every 3
                %pause(0.2);
            end
        end
    end
end