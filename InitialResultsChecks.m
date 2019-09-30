function [resultsChecksTable,usefulData]=InitialResultsChecks()
    dbstop if error

    expectedNumFiles=302;
    minNumRuns=0;
    maxCoV=0.05;
    
%     fileData=dir('results/*.mat');
%     numberOfFiles=size(fileData,1);
%     numFilesCorrect=expectedNumFiles==numberOfFiles;

    [filename,isPresent,isUpdated,numErrors,observedDifference,isObservedConverged,estimatedObjectiveCoV,isEstimatedConverged]=InitialiseVariables(expectedNumFiles);

    for batchIndex=1:expectedNumFiles
        [wInitS,wUltStrength,wPdStrain,wYStrength,wUltStrain] = findWeights(batchIndex);

        filename{batchIndex}=['results/MOBO/additive/pristine/MOBO_additive_pristine_results_' num2str(wInitS) '_' num2str(wUltStrain) '_' num2str(wPdStrain) '_' num2str(wUltStrength) '_' num2str(wYStrength) '.mat'];

        try
            load(filename{batchIndex})
        catch
            continue
        end

        isPresent(batchIndex)=true;
        isUpdated(batchIndex)=CheckBatchIndexUpdated(filename{batchIndex});
        [numErrors(batchIndex),observedDifference(batchIndex),isObservedConverged(batchIndex),estimatedObjectiveCoV(batchIndex),isEstimatedConverged(batchIndex)]=CheckResults(results,minNumRuns,maxCoV);

    end
    
    [resultsChecksTable,usefulData]=SummariseResults(expectedNumFiles,filename,isPresent,isUpdated,numErrors,observedDifference,isObservedConverged,estimatedObjectiveCoV,isEstimatedConverged);
    
    UpdateFileDataRecord()
end

function [filename,isPresent,isUpdated,numErrors,observedDifference,isObservedConverged,estimatedObjectiveCoV,isEstimatedConverged]=InitialiseVariables(expectedNumFiles)
    % intialises data at start to save time
    filename=cell(expectedNumFiles,1);
    isPresent=false(expectedNumFiles,1);
    isUpdated=false(expectedNumFiles,1);
    numErrors=zeros(expectedNumFiles,1);
    observedDifference=NaN(expectedNumFiles,1);
    isObservedConverged=false(expectedNumFiles,1);
    estimatedObjectiveCoV=NaN(expectedNumFiles,1);
    isEstimatedConverged=false(expectedNumFiles,1);
end

function [wInitS,wUltStrength,wPdStrain,wYStrength,wUltStrain] = findWeights(batchIndex)
    % finds the weights associated with each batch index (run number for HPC)
    wInitS=(batchIndex>=1 && batchIndex<=101).*(1-(batchIndex-1)./100);
    wUltStrength=(batchIndex>=1 && batchIndex<=101).*((batchIndex-1)./100)+(batchIndex>=102 && batchIndex<=201).*(2-(batchIndex-1)./100);
    wPdStrain=(batchIndex>101 && batchIndex<=201).*((batchIndex-1)./100-1);
    wYStrength=(batchIndex>=202 && batchIndex<=302).*(3-(batchIndex-2)./100);
    wUltStrain=(batchIndex>=202 && batchIndex<=302).*((batchIndex-2)./100-2);
    
    % scales these values to match filenames
    wInitS=wInitS*100;
    wUltStrength=wUltStrength*100;
    wPdStrain=wPdStrain*100;
    wYStrength=wYStrength*100;
    wUltStrain=wUltStrain*100;
end

function isUpdated=CheckBatchIndexUpdated(filepath)
    isUpdated=false;
    try 
        load('fileDataRecord.mat');
    catch
        % if file not found then then all data is new - could speed up
        % here!
        disp('no previous record of any file data')
        isUpdated=true;
        return
    end
    
    fileData=dir(filepath); 
    % remove results/ from filepath
    [~,filename]=strtok(filepath,'/');
    filename=filename(2:end);
    
    % finds the index where the filenames match
    dataRecordIndex=find(arrayfun(@(n) strcmp(fileDataRecord(n).name, filename), 1:numel(fileDataRecord)));
    if isempty(dataRecordIndex)==1
        disp(['no previous record of data for: ' filename])
        isUpdated=true;
        return
    end
    
    % if the date of the new file is later than the file data record, then
    % the data has been updated
    if datetime(fileData.date)>datetime(fileDataRecord(dataRecordIndex).date)
        isUpdated=true;
    end
    
end

function [numErrors,observedDifference,isObservedConverged,estimatedObjectiveCoV,isEstimatedConverged]=CheckResults(results,minNumRuns,maxCoV)
    % finds value of first evaluation of the objective function using the
    % evaluation time trace
    indexOfRunStart=find(isnan(results.ObjectiveEvaluationTimeTrace)==1,1,'last')+1;
    if isempty(indexOfRunStart)==1
        indexOfRunStart=1;
    end
    numRuns=results.NumObjectiveEvaluations-indexOfRunStart+1;
    
    % error trace gives 1 when objective function returns an error.
    numErrors=sum(results.ErrorTrace==1);
    
    % Observed minimum trace is thought ot have converged if no difference
    % in observed results
    observedDifference=results.ObjectiveMinimumTrace(end)-results.ObjectiveMinimumTrace(indexOfRunStart);
    isObservedConverged=(observedDifference==0) && (numRuns>=minNumRuns);
    
    % estimated objective is assumed to have convegred when the CoV of
    % estimated objected is less than 5%
    estimatedObjectiveCoV=std(results.EstimatedObjectiveMinimumTrace,'omitnan')/mean(results.EstimatedObjectiveMinimumTrace,'omitnan');
    isEstimatedConverged=abs(estimatedObjectiveCoV)<=maxCoV && (numRuns>=minNumRuns);
    
end

function [resultsChecksTable,usefulData]=SummariseResults(expectedNumFiles,filename,isPresent,isUpdated,numErrors,observedDifference,isObservedConverged,estimatedObjectiveCoV,isEstimatedConverged)
    resultsChecksTable=table(filename,isPresent,isUpdated,numErrors,observedDifference,isObservedConverged,estimatedObjectiveCoV,isEstimatedConverged);
    totalMissingResults=expectedNumFiles-sum(isPresent);
    totalMissingUpdates=expectedNumFiles-sum(isUpdated);
    totalErroneousRuns=sum(numErrors>0);
    totalObservedNotConverged=expectedNumFiles-sum(isObservedConverged);
    totalEstimatedObjectiveNotConverged=expectedNumFiles-sum(isEstimatedConverged);
    summaryData=struct('totalMissingResults',totalMissingResults,'totalMissingUpdates',totalMissingUpdates,'totalErroneousRuns',totalErroneousRuns,'totalObservedNotConverged',totalObservedNotConverged,'totalEstimatedObjectiveNotConverged',totalEstimatedObjectiveNotConverged);

    disp(['Total missing results: ' num2str(totalMissingResults)]);
    disp(['Total missing updates: ' num2str(totalMissingUpdates)]);
    disp(['Total erroneous runs: ' num2str(totalErroneousRuns)]);
    disp(['Total observed not converged: ' num2str(totalObservedNotConverged)]);
    disp(['Total estimated objective not converged: ' num2str(totalEstimatedObjectiveNotConverged)]);
    
    batchIndexMissingResults=find(isPresent==false);
    batchIndexMissingUpdates=find(isUpdated==false);
    batchIndexErroneousRuns=find(numErrors>0);
    batchIndexObservedNotConverged=find(isObservedConverged==false);
    batchIndexEstimatedObjectiveNotConverged=find(isEstimatedConverged==false);
    batchIndexData=struct('batchIndexMissingResults',batchIndexMissingResults,'batchIndexMissingUpdates',batchIndexMissingUpdates,'batchIndexErroneousRuns',batchIndexErroneousRuns,'batchIndexObservedNotConverged',batchIndexObservedNotConverged,'batchIndexEstimatedObjectiveNotConverged',batchIndexEstimatedObjectiveNotConverged);
    
    filenamesMissingResults=filename(isPresent==false);
    filenamesMissingUpdates=filename(isUpdated==false);
    filenamesErroneousRuns=filename(numErrors>0);
    filenamesObservedNotConverged=filename(isObservedConverged==false);
    filenamesEstimatedObjectiveNotConverged=filename(isEstimatedConverged==false);
    filenameData=struct('filenamesMissingResults',{filenamesMissingResults},'filenamesMissingUpdates',{filenamesMissingUpdates},'filenamesErroneousRuns',{filenamesErroneousRuns},'filenamesObservedNotConverged',{filenamesObservedNotConverged},'filenamesEstimatedObjectiveNotConverged',{filenamesEstimatedObjectiveNotConverged});

    usefulData=struct('summaryData',summaryData,'batchIndexData',batchIndexData,'filenameData',filenameData);
end

function UpdateFileDataRecord()
    % update the file data record for next time
    fileDataRecord=dir('results/MOBO/additive/pristine/*.mat');
    save('matlab_scripts/fileDataRecord','fileDataRecord','-v7.3');
end