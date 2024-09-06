module App
using Main.StatisticAnalysis
using GenieFramework
using GenieCache
using PlotlyBase
# using PlotlyJS
# using Plots
using Statistics
using KernelDensity
using Interpolations
using DataFrames
using Colors
using CSV
using BSON
using StippleDownloads
# Plots.plotly()
@genietools



@app begin
    @out fileIndex=1 # store which experiment we are on
    @out dataType="stained" # first will be stained and then unstained 
    @in xlabel ="" # labels for the axes
    @in ylabel ="" # labels for the axes
    @out p="upload first sample/cell line" # text used on the upload box so we can change our instructions to the user
    const FILE_PATH = joinpath("public", "uploads") # This is where user uploads go
    mkpath(FILE_PATH)
    @out upfiles = [] # I'm not sure we use this TODO: check deleting this. 
    @out thisLayout = PlotlyBase.Layout(title="Fingerprint") # The layout of the main fingerprint plot
    # A secondary plot only used to show a legend
    @out thisLayout2 = PlotlyBase.Layout(title="",yaxis_visible=false,xaxis_visible=false,showlegend = true,height=200,plot_bgcolor = "white")
    @private data = DataFrame() # Data is stored in a dataframe once uploaded
    @in selected_file = ""# This stores the last file uploaded,
    @in selected_column1 = ""# This stores the column selected to plot
    @in selected_column2 = ""# This stores the other column we're plotting
    @out columns = [] # store all columns
    
    # Controls if we want to highlight the upload box, or one of the two axes naming text boxes
    @out uploadBoolean=true
    @out highlightCol1=false
    @out highlightCol2=false
    @out highlightLabel1=true
    @out highlightLabel2=false
    
    # Controls whether the last selections is visible, we hide it in some scenarios
    @out hideSelections1=false
    @out hideSelections2=false
    
    # Store which column was selected as the one to plot
    @out columnSelectText1 = ""
    @out columnSelectText2 = ""
    @out instructionText = "" # Present some dynamic instruction text on the main screen.
    @out columnDict = Dict() # Dictionary of columns
    
    # The name of the last stained experiment needs to be stored, this is important because the last file uploaded
    # may be unstained but we want to identify the name of an experiment by the stained experiment file name
    @out lastStained=""
    
    @out stainedUnstainedDict = Dict() # A dictionary to get matching unstained experiment from the name of the stained one
    
    @out allDatasets = []# A list of all datasets (stained file names), to let us pick which one needs editing in the edit menu
    
    # Once data has been uploaded and columne selected we want to disable the column selection and naming boxes
    @out namingBoolean = false
    @out columnBoolean = true

    @out fileUploadLabel = "" # dynamic text for the upload box
    @in dimensions = "2D" # default to 2D but change if needed
    @in normalizeHeight = "No" # default to non-normalised but change if needed
    @out toggleLabel = "Dimensions: 2D" # dynamic label text to tell the user which type of plot we're making
    @out normalizeLabel = "Normalise: No"
    
    @in line_width= 4 # width of the lines contour lines
    @in num_contours=5 # how many contours (note this is total, may not be this may visible for all lines)
    @in mid_point=0.1 # Controls where the fade to white from colour reaches its mid point
    
    @out trace = []
    @out trace2 = []
    # @out trace = []
    @in enable2DControls=false # some controls are only available in 2D
    @in colors=["#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","gray"] #default colour scheme
    @in Range_x = RangeData(30:60) #default range for both axes
    @in Range_y = RangeData(30:60)
    @in customHex = "#CC79A7" #placeholder text for color textbox
    @in lineToEdit = "select experiment" # Placeholder text for fingerprint editing selection box
    @in colorBoolean = true # color editing is only enabled when a line is selected
    @in linecolorDict=Dict() # link between experiment numbers and their colours
    @in delete_line = false # line deleting is only available when a line is selected
    @out toggleCentreLabel = "Add centre dot: No" # text that tells user the status of middle dots
    @in centreBoolean = "No" 
    @in save_all = false # triggered when saveall button is pressed
    @in purge_cache = false # triggered when purge_cache (start again) button is pressed
    @in downloadButtonBoolean = true # download BSON button only valid sometimes
    @out completeDataDictStained=Dict() # all the stained data is stored here so we don't need to read CSV files multiple times
    @out completeDataDictUnstained=Dict()# all the unstained data is stored here so we don't need to read CSV files multiple times
    @in updatingManyVars = false
    @in tab_selected = "Upload" # this is used to store which tab the user is currently in, we start in upload.

    # we need a function to update the dataframe that holds all the data. This is where the CSV will happen.
    # if a user has loaded a session from a BSON file the CSV file uploads won't exist anymore so this should only be called
    # when an upload has happened.
    function updateDFwithLatestCSVs(lastStained,FILE_PATH,completeDataDictStained,completeDataDictUnstained,stainedUnstainedDict)
        # look up the name of unstained and stained experiments
        thisStained=lastStained
        thisUnstained=stainedUnstainedDict[lastStained]
        #read both CSV files
        dataStained = CSV.read(joinpath(FILE_PATH, thisStained), DataFrame)
        dataUnstained = CSV.read(joinpath(FILE_PATH, thisUnstained), DataFrame)
        #add the dataframes created from these CSV files to the main data dictionaries
        completeDataDictStained[thisStained]=dataStained
        completeDataDictUnstained[thisStained]=dataUnstained
        # update the main data dictionaries
        return completeDataDictStained,completeDataDictUnstained
    end

    # when any changes happen we want to update the plot, this means going through the data dictionary and replotting every line with
    # the latest parameters. This function returns two sets of layouts and traces (one for the main plot and one for the legend)
    # every paramter needed to make the plot is passed in.
    function updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)

        # we need to merge multiple experiments into a single dataframe
        # this lets us calculate Z-scores
        mergedDF=DataFrame([String[],Float64[],Float64[]], ["experiment", xlabel,ylabel])

        # We want to make a fingerprint for every stained experiment
        allStained=collect(keys(stainedUnstainedDict))

        # for each experiment
        for i in 1:length(allStained)
            # get the name of the stained and unstained experiemnt
            thisStained=allStained[i]
            thisUnstained=stainedUnstainedDict[allStained[i]]
            # look up the data for these experiments in the main data dictionary
            dataStained=completeDataDictStained[allStained[i]]
            dataUnstained=completeDataDictUnstained[allStained[i]]
            
            # We need to put all events into an array to take medians etc
            numberOfEvents=size(dataStained,1)
            ArrayToPopulate=Array{Any}(undef, numberOfEvents, 3)
            # store the x and y values for the stained and unstained
            # TODO: ensure all this variable swapping around is necessary...
            (column1Stained,column2Stained)=columnDict[thisStained]
            (column1UnStained,column2UnStained)=columnDict[thisUnstained]
            eachMeasurement1=dataStained[:,column1Stained]
            eachMeasurement2=dataStained[:,column2Stained]
            eachMeasurementUnStained1=dataUnstained[:,column1UnStained]
            eachMeasurementUnstained2=dataUnstained[:,column2UnStained]

            # take the median of both unstained data columns
            thisUnstainedMFI1=median(eachMeasurementUnStained1)
            thisUnstainedMFI2=median(eachMeasurementUnstained2)
            # substract the median of the unstained from each stained data value.
            ToAdd1=eachMeasurement1.-thisUnstainedMFI1
            ToAdd2=eachMeasurement2.-thisUnstainedMFI2
            # we need to store each value in one large array, with the name of the experiment in the first column
            # this lets us take z-score of the dataset as whole, while keeping track of which experiment is which
            ArrayToPopulate[:,1].=thisStained
            ArrayToPopulate[:,2].=Float64.(ToAdd1)
            ArrayToPopulate[:,3].=Float64.(ToAdd2)
            ArrayToPopulate=DataFrame(ArrayToPopulate,names(mergedDF))
            # zeros mess things up so filter out anything where the sum of both data columns isn't two.
            # TODO: this is not the best way to do this.
            ArrayToPopulate=filter(row -> (sum([row[2],row[3]])>2),  ArrayToPopulate)
            # add this experiment onto the large merged dataframe of all experiments
            append!(mergedDF,ArrayToPopulate)
        end
        
        # we save the output on the server, at some point I wanted to let people download this but we don't currently
        CSV.write("testOutput.csv", mergedDF)
        # calculate the z-score
        # subtract the mean from each datapoint and then divide by the standard deviation
        meanCol2=mean(mergedDF[:,2]) 
        stdCol2=std(mergedDF[:,2])
        meanCol3=mean(mergedDF[:,3]) 
        stdCol3=std(mergedDF[:,3])
        mergedDF[:,2].= (mergedDF[:,2] .- meanCol2) ./ stdCol2
        mergedDF[:,3].= (mergedDF[:,3] .- meanCol3) ./ stdCol3
        # write the z-scored data out as well (we don't use this)
        CSV.write("testOutputStd.csv", mergedDF)

        # we need to build up a list of traces, one for each fingerprint (and centre dot)
        trace=[]
        trace2=[]
        # we also need to define two layouts, one for the main plot and one for the legend
        # we do this so the legend doesn't block any data, whatever the data
        thisLayout=PlotlyBase.Layout()
        thisLayout2=PlotlyBase.Layout()

        # for every experiment
        for i in 1:length(allStained)
            # get the part of the big merged and standardised dataframe that represents this experiment
            thisDF=mergedDF[(mergedDF.experiment .== allStained[i]), :]
            # calculate the kernel density of the points in this fingerprint
            k = kde((thisDF[:,2], thisDF[:,3]),npoints=(100,100),boundary=((-5,5),(-5,5)))
            # if we want to we can normalise each KDE to max out at 1
            if normalizeHeight == "Yes"
              k.density=k.density./maximum(k.density)
            end
            
            # the way we plot 2D and 3D plots is totally different
            if dimensions == "3D"
                # we're working in 3D so we need a surface plot
                # TODO: let people change the opactity
                thisTrace = surface(z=k.density, line_width=line_width,opacity=0.3,showscale=false,colorscale=[[0, "white"], [1, linecolorDict[allStained[i]]]],contours_z=attr(
                    show=true,
                    usecolormap=false,
                    # colour by which experiment they belong to
                    color=linecolorDict[allStained[i]],
                    highlightcolor="limegreen",
                    # this potion puts the contours on the floor or ceiling
                    project_z=true,
                    line_width=line_width,
                    z=attr(show=true, start= 0.1, size= 0.1),
                    z_end=0.5
                ))
                
                # add this surface to the list of surfaces
                trace=vcat(trace,[thisTrace])
                thisLayout = PlotlyBase.Layout(title=xlabel*" against "*ylabel*" fingerprint",
                # set the range and labels appropriately.
                scene=attr(xaxis_title=xlabel,
                yaxis_title=ylabel, 
                zaxis_title="density"),
                scene_xaxis_range=[Range_x.range.start, Range_x.range.stop],
                scene_yaxis_range=[Range_y.range.start, Range_y.range.stop],
                scene_zaxis_range=[0.1, 1.0],
                showlegend = false,
                showscale = false,
                height=600
                )
            else
                # we are in two dimensions, so we need a contour plot
                # we will make the colour scale fade from white to the appropriate colour for this experiement.
                # the user can control the mid point of this fade to make it half way or at an extreme such that everyline is coloured
                thisTrace = contour(;z=k.density,line_width=line_width,colorscale=[[0, "white"],[mid_point,linecolorDict[allStained[i]]],[1, linecolorDict[allStained[i]]]], contours_coloring="lines",showscale=false,ncontours=num_contours)
                # add this contour to the list of lines to plot
                trace=vcat(trace,[thisTrace])
                
                # do we want to add a dot in the middle?
                if centreBoolean =="Yes"
                    # if so put this dot at the MFI of each variable
                    yCentre=median(thisDF[:,2])
                    xCentre=median(thisDF[:,3])
                    # because we did the density plot from -5 to 5 (a range of 10) with 100 points we need to start  with the centre at 50
                    # and scale it 10x higher than whatever value it has
                    # So if the MFI of the z-score is at 1, we need to plot it at (50 + 10*1) = 60
                    # if the MFI of the z-score is -1 we need to plot it at (50 + 10*-1) = 40
                    xCentre=50+(10*xCentre)
                    yCentre=50+(10*yCentre)
                    # we need to plot the point and add it to our list of traces
                    extraPoint=scatter(x=[xCentre,xCentre],y=[yCentre,yCentre],mode="markers",color=linecolorDict[allStained[i]],marker=attr(size=12,color=linecolorDict[allStained[i]]))
                    trace=vcat(trace,[extraPoint])
                end

                # now we have all traces for all experiments so lets just style the overall plot
                thisLayout = PlotlyBase.Layout(title=xlabel*" against "*ylabel*" fingerprint",xaxis_title=xlabel,yaxis_title=ylabel,
                scene=attr(xaxis_title=xlabel,
                yaxis_title=ylabel, 
                zaxis_title="density"),
                xaxis_range=[Range_x.range.start, Range_x.range.stop],
                yaxis_range=[Range_y.range.start, Range_y.range.stop],
                yaxis=attr(linecolor="black",gridcolor="lightgrey",ticks="outside"),
                xaxis=attr(linecolor="black",gridcolor="lightgrey",ticks="outside"),
                plot_bgcolor = "white",
                showscale = false,
                showlegend = false,                
                height=600)

            end
            # We also need to add an arbitrary trace and a legend to our second plot. 
            # this lets us have a legend we can export and arrange wherever we want in illustrator
            thisTrace2=scatter(x=[1,2],y=[1,1],size=(2,2),mode="lines",line_width=line_width,line=attr(color=linecolorDict[allStained[i]],width=line_width),name=allStained[i])
            trace2=vcat(trace2,[thisTrace2])
            thisLayout2=PlotlyBase.Layout(title="legend",yaxis_visible=false,xaxis_visible=false,showlegend = true,height=200,plot_bgcolor = "white")
        end
        # return the trace and layout for both plots which will update the plots dynamically.
        return trace, thisLayout,trace2,thisLayout2
    end

    @onchange dimensions begin    
        # when the user change the dimension slider
        if !updatingManyVars # make sure we're not loading a dataset currently
            # change the dimension variable and re-run update plot which will replot everything correctly
            toggleLabel= "Dimensions: " *dimensions
            trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
            enable2DControls=(dimensions=="3D")
        end
    end

    @onchange centreBoolean begin
        # when a user adds or removed the centre point
        if !updatingManyVars # make sure we're not loading a dataset currently
            # todo: make these actual booleans, this makes me sad!
            if centreBoolean == "Yes"
                toggleCentreLabel= "Add centre dot: Yes"
            else
                toggleCentreLabel= "Add centre dot: No"
            end
            # simply replot the plot now that the boolean has changed
            trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
        end
end

    @onchange lineToEdit begin    
        # the user has selected a line to edit the colour of or delete
        if !updatingManyVars # make sure we're not loading a dataset currently
            # allow the user to edit colours of plots now
            colorBoolean=false
        end
    end

    @onchange customHex begin
        # the user has input a custom colour
        if !updatingManyVars # make sure we're not loading a dataset currently
            # change the colour to match what the user chose and rerun the plotting code    
            linecolorDict[lineToEdit]=customHex
            trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
            
        end
    end

    @onbutton delete_line begin
        # the user has pressed "delete" on a line
        if !updatingManyVars # make sure we're not loading a dataset currently
            # simply remove thie line from the stained/Unstained lookup dicionary
            # this is used to loop through all experiments in the plotting code so the experiment
            # vanishes from the plot. Then update the plot of course!
            delete!(stainedUnstainedDict,lineToEdit)
            trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
        end
    end

    @event :purge_cache begin  
        # the user wants to start again
        # a refresh works in dev but not prod
        # Here we set all variables back to their default value

        # Tell all the other reactive variable sthat that we're about to make many changes
        # this stops the plot trying to regenerate throughout this funciton
        updatingManyVars = true

        # set all variables to their initial value (there must be a better way)
        fileIndex=1
        dataType="stained"
        xlabel =""
        ylabel =""
        p="upload first sample/cell line"
        upfiles = []
        thisLayout = PlotlyBase.Layout(title="Fingerprint")
        thisLayout2 = PlotlyBase.Layout(title="",yaxis_visible=false,xaxis_visible=false,showlegend = true,height=200,plot_bgcolor = "white")
        data = DataFrame()
        selected_file = ""
        selected_column1 = ""
        selected_column2 = ""
        columns = []
        uploadBoolean=true
        highlightCol1=false
        highlightCol2=false
        highlightLabel1=true
        highlightLabel2=false
        hideSelections1=false
        hideSelections2=false
        columnSelectText1 = ""
        columnSelectText2 = ""
        instructionText = ""
        columnDict = Dict()
        lastStained=""
        stainedUnstainedDict = Dict()
        allDatasets = []
        namingBoolean = false
        columnBoolean = true
        fileUploadLabel = ""
        dimensions = "2D"
        toggleLabel = "Dimensions: 2D"
        normalizeLabel = "Normalise: No"
        normalizeHeight = "No"
        line_width= 4
        num_contours=5
        mid_point=0.1
        trace = []
        trace2 = []
        enable2DControls=false
        colors=["#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","gray"]
        Range_x = RangeData(30:60)
        Range_y = RangeData(30:60)
        customHex = "#CC79A7"
        lineToEdit = "select experiment"
        colorBoolean = true
        linecolorDict=Dict()
        delete_line = false
        toggleCentreLabel = "Add centre dot: No"
        centreBoolean = "No"
        save_all = false
        purge_cache = false
        downloadButtonBoolean = true
        completeDataDictStained=Dict()
        completeDataDictUnstained=Dict()
        updatingManyVars = false

        # purge Genie's cache (this may not be necessary)
        GenieCache.purge()
    end
    
    @event save_all begin
        # user has pressed the save workspace button
        # we will make copies of all variables, because you can't save reactive variables
        # to a BSON file. Then we will save these copies to a BSON file. When we load we
        # do the reverse

        # reassure the user the download is coming because this can take a second or two
        notify(__model__,"preparing download")

        # get ready to write out the BSON
        io = IOBuffer()

        # we have to create copies of the variables to be able to save them with BSON :(
        data_fileIndex = fileIndex
        data_dataType = dataType
        data_xlabel=xlabel
        data_ylabel=ylabel
        data_p=p
        data_upfiles = upfiles
        data_thisLayout=thisLayout
        data_thisLayout2=thisLayout2
        data_data=data
        data_selected_file=selected_file
        data_selected_column1=selected_column1
        data_selected_column2=selected_column2
        data_columns=columns
        data_uploadBoolean=uploadBoolean
        data_highlightCol1=highlightCol1
        data_highlightCol2=highlightCol2
        data_highlightLabel1=highlightLabel1
        data_highlightLabel2=highlightLabel2
        data_hideSelections1=hideSelections1
        data_hideSelections2=hideSelections2
        data_columnSelectText1=columnSelectText1
        data_columnSelectText2=columnSelectText2
        data_instructionText=instructionText
        data_columnDict=columnDict
        data_lastStained=lastStained
        data_stainedUnstainedDict=stainedUnstainedDict
        data_allDatasets=allDatasets
        data_namingBoolean=namingBoolean
        data_columnBoolean=columnBoolean
        data_fileUploadLabel=fileUploadLabel
        data_dimensions=dimensions
        data_toggleLabel=toggleLabel
        data_normalizeLabel=normalizeLabel
        data_normalizeHeight=normalizeHeight
        data_line_width=line_width
        data_num_contours=num_contours
        data_mid_point=mid_point
        data_trace=trace
        data_trace2=trace2
        data_enable2DControls=enable2DControls
        data_colors=colors
        data_Range_x=Range_x
        data_Range_y=Range_y
        data_customHex=customHex
        data_lineToEdit=lineToEdit
        data_colorBoolean=colorBoolean
        data_linecolorDict=linecolorDict
        data_delete_line=delete_line
        data_toggleCentreLabel=toggleCentreLabel
        data_centreBoolean=centreBoolean
        data_save_all=save_all
        data_completeDataDictStained=completeDataDictStained
        data_completeDataDictUnstained=completeDataDictUnstained

        # save all the variables to a BSON file
        BSON.@save io data_fileIndex data_dataType data_xlabel data_ylabel data_p data_upfiles 	data_thisLayout	data_thisLayout2	data_data	data_selected_file	data_selected_column1	data_selected_column2	data_columns	data_uploadBoolean	data_highlightCol1	data_highlightCol2	data_highlightLabel1	data_highlightLabel2	data_hideSelections1	data_hideSelections2	data_columnSelectText1	data_columnSelectText2	data_instructionText	data_columnDict	data_lastStained	data_stainedUnstainedDict	data_allDatasets	data_namingBoolean	data_columnBoolean	data_fileUploadLabel	data_dimensions	data_toggleLabel	data_normalizeLabel	data_normalizeHeight	data_line_width	data_num_contours	data_mid_point	data_trace	data_trace2	data_enable2DControls	data_colors	data_Range_x	data_Range_y	data_customHex	data_lineToEdit	data_colorBoolean	data_linecolorDict	data_delete_line	data_toggleCentreLabel	data_centreBoolean	data_save_all data_completeDataDictStained data_completeDataDictUnstained
        seekstart(io)
        # trigger a download for the user.
        try
            download_binary(__model__, take!(io), "data.bson")
        catch ex
            println(ex)
            notify(__model__,"error generating download"*ex)
        end
    end
    

    @onchange line_width, num_contours, mid_point, Range_x, Range_y begin    
        # There are many parameters that simply need replotting when they're changed so let's do just that.
        if !updatingManyVars # make sure we're not loading a dataset currently
            trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
        end
    end
    
    @onchange normalizeHeight begin
        # user toggled normalise. Here we need to change some text about the status of the toggle and replot
        if !updatingManyVars # make sure we're not loading a dataset currently
            normalizeLabel= " Normalise: " * normalizeHeight
                trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
        end
    end
    @onchange selected_column1 begin
        # the user has selected the first column in their dataset
        if !updatingManyVars # make sure we're not loading a dataset currently
            # highlight tht they can now select a second column and deemphasise the first column
            highlightCol1 = false
            highlightCol2 = true
            hideSelections1=false
        end
        
    end
    

    @onchange selected_column2 begin
        # the user has selected the second colum in their dataset and is now ready to either:
        # plot the next dataset if this was an unstained control
        # or upload the unstained control if this was a stained experiment.
        if !updatingManyVars # make sure we're not loading a dataset currently
            # change what is focused now
            hideSelections2=false
            highlightCol2 = false
            columnBoolean=true
            # add these columns to the column dictionary
            columnDict[selected_file]=(selected_column1,selected_column2)
            
            # Check whether this was stained or unstained
            if dataType == "stained" 
                # if it was stained it is time to ask for a matching unstained
                dataType="unstained"
                instructionText = "Return to file upload and upload the control experiment matching this sample."
                uploadBoolean = false
                lastStained=selected_file
                fileUploadLabel="Upload control CSV file for this dataset:"
            else 
                # if this was unstained it's time to ask the user whether they want to upload another dataset
                dataType="stained"
                notify(__model__,"Now upload next sample if there is one")
                # we now have one more complete dataset
                fileIndex=fileIndex+1
                instructionText = "Return to file upload and upload the next sample if there is one."
                fileUploadLabel="Upload fully stained CSV file for the next dataset:"
                uploadBoolean = false
                # save which file is the stained file and which is the unstained pairing in the dictionary
                stainedUnstainedDict[lastStained]=selected_file
                allDatasets=collect(keys(stainedUnstainedDict))
                # at this point we have all the data we need to make or update a fingerprint

                # if there are enough colours in the colour list asign the next colour to this line, otherwise black. 
                try
                    linecolorDict[lastStained]=colors[length(linecolorDict)+1]
                catch e
                    notify(__model__,"We've run out of colours, assuming black.")
                    linecolorDict[lastStained]="#000000"
                end
                
                # we can now update our big data dictionary of all stained and unstained data with these new CSV files
                completeDataDictStained,completeDataDictUnstained=updateDFwithLatestCSVs(lastStained,FILE_PATH,completeDataDictStained,completeDataDictUnstained,stainedUnstainedDict)
                # update the plot now that the new data is in the dictionary and it should appear
                trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)                
            end
        end
    end
    
    @onchange selected_file begin
        # user has uploaded a new file
        if !updatingManyVars # make sure we're not loading a dataset currently
            # read the file and get its columns so the user can choose which ones are needed here
            data = CSV.read(joinpath(FILE_PATH, selected_file), DataFrame)
            columns = names(data)
            # allow the user to select columns but not do anything else at this point.
            uploadBoolean=true
            namingBoolean = true
            columnBoolean = false
            hideSelections1=true
            hideSelections2=true
        end
    end

    @onchange xlabel,ylabel,fileuploads begin
        # The user has changed labels or we have an experiment name now. 
        # we don't need a total replot just a refresh of the layout
        if !updatingManyVars # make sure we're not loading a dataset currently
            layout = PlotlyBase.Layout(title=xlabel*" against "*ylabel*"",height=800)
        end
    end

    @onchange xlabel begin
        # the user has given the x axis a name
        if !updatingManyVars # make sure we're not loading a dataset currently
            # We highlight naming the next axis and change the text on the column selection box
            columnSelectText1="select column for:"*xlabel
            highlightLabel1=false
            highlightLabel2=true
        end
    end
    @onchange ylabel begin
        # the user has given the x axis a name
        if !updatingManyVars # make sure we're not loading a dataset currently
            # We highlight naming the next axis and change the text on the column selection box
            # As this is the second one we now want them to upload a file
            columnSelectText2="select column for:"*ylabel
            uploadBoolean = false
            fileUploadLabel = "Upload fully stained CSV file for this dataset:"
            highlightLabel1=false
            highlightLabel2=false
        end
    end

    @onchange fileuploads begin
        # a file upload has happened
        if ! isempty(fileuploads) && !updatingManyVars
            @info "File was uploaded: " fileuploads
            filename = fileuploads["name"]
            # the file will either be a BSON (reloading an entire workspace)
            # or a CSV file, we need to process them similarly but different
            if last(filename,4)=="bson" || last(filename,4)=="BSON"
                # if it's a BSON file we just need to move it to the uploads folder/
                try
                    isdir(FILE_PATH) || mkpath(FILE_PATH)
                    mv(fileuploads["path"], joinpath(FILE_PATH, filename), force=true)
                catch e
                    println("ERROR! PATH: "*fileuploads["path"])
                    @error "Error processing file: $e"
                    notify(__model__,"Error processing file: $(fileuploads["name"])")
                end
                thisBSON=joinpath(FILE_PATH, filename)

                # we want to disable most responses to reactive variables at the point because we're about to change many variables
                updatingManyVars = true
                # load every variable in the BSON
                BSON.@load thisBSON data_fileIndex data_dataType data_xlabel data_ylabel data_p data_upfiles data_thisLayout	data_thisLayout2	data_data	data_selected_file	data_selected_column1	data_selected_column2	data_columns	data_uploadBoolean	data_highlightCol1	data_highlightCol2	data_highlightLabel1	data_highlightLabel2	data_hideSelections1	data_hideSelections2	data_columnSelectText1	data_columnSelectText2	data_instructionText	data_columnDict	data_lastStained	data_stainedUnstainedDict	data_allDatasets	data_namingBoolean	data_columnBoolean	data_fileUploadLabel	data_dimensions	data_toggleLabel	data_normalizeLabel	data_normalizeHeight	data_line_width	data_num_contours	data_mid_point	data_trace	data_trace2	data_enable2DControls	data_colors	data_Range_x	data_Range_y	data_customHex	data_lineToEdit	data_colorBoolean	data_linecolorDict	data_delete_line	data_toggleCentreLabel	data_centreBoolean	data_save_all data_completeDataDictStained data_completeDataDictUnstained

                # remember we changed the names to make sure we were working with copies of reactive variables?
                # now we undo that.
                completeDataDictStained = data_completeDataDictStained
                completeDataDictUnstained= data_completeDataDictUnstained
                upfiles=data_upfiles 
                thisLayout=data_thisLayout
                thisLayout2=data_thisLayout2
                p=data_p
                data=data_data
                fileIndex=data_fileIndex 
                dataType=data_dataType 
                xlabel=data_xlabel
                ylabel=data_ylabel
                selected_file=data_selected_file
                selected_column1=data_selected_column1
                selected_column2=data_selected_column2
                columns=data_columns
                uploadBoolean=data_uploadBoolean
                highlightCol1=data_highlightCol1
                highlightCol2=data_highlightCol2
                highlightLabel1=data_highlightLabel1
                highlightLabel2=data_highlightLabel2
                hideSelections1=data_hideSelections1
                hideSelections2=data_hideSelections2
                columnSelectText1=data_columnSelectText1
                columnSelectText2=data_columnSelectText2
                instructionText=data_instructionText
                columnDict=data_columnDict
                lastStained=data_lastStained
                stainedUnstainedDict=data_stainedUnstainedDict
                allDatasets=data_allDatasets
                namingBoolean=data_namingBoolean
                columnBoolean=data_columnBoolean
                fileUploadLabel=data_fileUploadLabel
                dimensions=data_dimensions
                toggleLabel=data_toggleLabel
                normalizeLabel=data_normalizeLabel
                normalizeHeight=data_normalizeHeight
                line_width=data_line_width
                num_contours=data_num_contours
                mid_point=data_mid_point
                trace=data_trace
                trace2=data_trace2
                enable2DControls=data_enable2DControls
                colors=data_colors
                Range_x=data_Range_x
                Range_y=data_Range_y
                customHex=data_customHex
                lineToEdit=data_lineToEdit
                colorBoolean=data_colorBoolean
                linecolorDict=data_linecolorDict
                delete_line=data_delete_line
                toggleCentreLabel=data_toggleCentreLabel
                centreBoolean=data_centreBoolean
                save_all=data_save_all
                # okay we are done and we want our variables to be reactive again
                updatingManyVars = false

                # update the plot with all the loaded variables!
                trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean,completeDataDictStained,completeDataDictUnstained)
            else
                # this must be a CSV file
                try
                    isdir(FILE_PATH) || mkpath(FILE_PATH)
                    # move it to the upload folder but give it a name that includes whether it's stained or unstained and which experiment it is
                    mv(fileuploads["path"], joinpath(FILE_PATH, string(fileIndex)*"_"*dataType*"_"*filename), force=true)
                catch e
                    @error "Error processing file: $e"
                    notify(__model__,"Error processing file: $(fileuploads["name"])")
                    println(e)
                end
                newFileName=string(fileIndex)*"_"*dataType*"_"*filename
                # save that this new file exists, which will trigger the appropriate event
                selected_file=newFileName
                
            end
            # empty this dictionary as the event is over
            fileuploads = Dict{AbstractString,AbstractString}()
        end
        # a list of all uploaded files TODO: check we need this?
        upfiles = readdir(FILE_PATH)
    end

    
    @event uploaded begin
        # an upload has happened, we now what the user to focus on selecting their important columns in the data
        @info "uploaded"
        notify(__model__, "File was uploaded. Now select column.")
        columnBoolean = false
        highlightCol1 = true
        fileUploadLabel=""
        uploadBoolean=true
        updatingManyVars = true
        selected_column1 =""
        selected_column2 =""
        updatingManyVars = false
    end
    @event rejected begin
        # this happened automatically if the user uploads a non CSV/BSON file or tries to upload the same file twice.
        @info "rejected"
        notify(__model__, "Please upload a valid file")
    end

    @event uploadedBSON begin
        # the upload was a BSON, so we can notify the user
        @info "uploaded"
        notify(__model__, "previous data loaded")
        
    end
    @event rejectedBSON begin
        # if the BSOn was rejected let's tell the user
        @info "rejected"
        notify(__model__, "Please upload a valid file")
    end
    
end



function ui()
    # ther overall layout is two columns, one for uploading and settings, the other for the plots
    row([
        cell(class="col-md-5", [
            # Have three tabs to choose from
            tabgroup(
                :tab_selected,
                inlinelabel = true,
                class = "bg-primary text-white shadow-2",
                [
                    tab(name = "Upload", icon = "upload_file", label = "Upload"),
                    tab(name = "Edit", icon = "edit", label = "Edit"),
                    tab(name = "Tune", icon = "tune", label = "Tune"),
                ],
            ),
            tabpanels(
                :tab_selected,
                animated = true,
                var"transition-prev" = "scale",
                var"transition-next" = "scale",
                [
                    tabpanel(name = "Upload", [ 
                        # The upload tab
                        p("Upload your data, one experiment at a time."),
                        textfield("First protein name:", :xlabel,placeholder="Bcl2",disable=:namingBoolean,standout=:highlightLabel1),
                        textfield("Second protein name:", :ylabel,placeholder="Bclxl",disable=:namingBoolean,standout=:highlightLabel2),
                        uploader( multiple = false,
                            accept = ".csv",
                            maxfilesize = 1024*1024*150, # bytes
                            maxfiles = 30,
                            autoupload = true,
                            hideuploadbtn = true,
                            label = :fileUploadLabel,
                            nothumbnails = true,
                            style="max-width: 95%; width: 95%; margin: 0 auto;",
                            disable=:uploadBoolean,
                            @on("rejected", :rejected),
                            @on("uploaded", :uploaded)
                            ),
                        p("{{columnSelectText1}}"),
                        Stipple.select(:selected_column1; options=:columns,disable=:columnBoolean,standout=:highlightCol1,usechips=true,clearable=false,hideselected=:hideSelections1),
                        p("{{columnSelectText2}}"),
                        Stipple.select(:selected_column2; options=:columns,disable=:columnBoolean,standout=:highlightCol2,usechips=true,clearable=false,hideselected=:hideSelections2),
                        p("{{instructionText}}")
                    ]),
                    tabpanel(name = "Edit", [
                        # The edit panel
                        list(
                            bordered = true,
                            separator = true,
                            [
                            item_section(
                            [
                                item_label("Change color or delete an existing fingerprint."),
                                item(
                                    Stipple.select(:lineToEdit; options=:allDatasets,disable=false,usechips=true,clearable=true)),
                                item(
                                    textfield("Color (as hex code):", :customHex,placeholder="#CC79A7",disable=:colorBoolean)),
                                item(
                                    btn("Delete", icon = "delete", color = "red", class = "q-mr-sm",@click(:delete_line))),
                            ]), 
                            separator(spaced = true),
                            item_section(
                            [
                                item_label("Save entire workspace."),
                                item([
                                    btn("Save", icon = "cloud_download", color = "green", class = "q-mr-sm",@on(:click,:save_all,:addclient))]
                                ),
                                
                            ]),
                            separator(spaced = true),
                            item_section(
                            [
                                item_label("Load entire workspace."),
                                item(
                                    uploader( multiple = false,
                                    accept = ".bson",
                                    maxfilesize = 1024*1024*150, # bytes
                                    maxfiles = 1,
                                    autoupload = true,
                                    hideuploadbtn = true,
                                    label = "Upload a BSON file saved previously",
                                    nothumbnails = true,
                                    style="max-width: 95%; width: 95%; margin: 0 auto;",
                                    @on("rejected", :rejectedBSON),
                                    @on("uploaded", :uploadedBSON)
                                    ),
                                    
                                ),
                                
                            ]),
                            separator(spaced = true),
                            item_section(
                            [
                                item_label("Start again"),
                                item([
                                    # btn("Generate Downloadable File", icon = "convert_to_text", color = "blue", class = "q-mr-sm",@click(:save_all)),
                                    btn("Purge", icon = "cached", color = "grey", class = "q-mr-sm",@on(:click,:purge_cache))]
                                ),
                                
                            ])
                            ]
                        )]
                        ),
                    

                
                    tabpanel(name = "Tune", [
                        # The tune (visual settings) panel
                        list(
                            bordered = true,
                            separator = true,
                            [
                            item(
                            toggle(
                                :toggleLabel,
                                :dimensions,
                                leftlabel = true,
                                # label= :toggleLabel,
                                color = "blue",
                                val= 2,
                                var"true-value" = "3D",
                                var"false-value" = "2D",
                            )),
                            item(
                            toggle(
                                :normalizeLabel,
                                :normalizeHeight,
                                leftlabel = true,
                                # label= :normalizeLabel,
                                color = "blue",
                                # val= 2,
                                var"true-value" = "Yes",
                                var"false-value" = "No",
                            )),
                            item([
                            # p("Line Width:"),
                            itemsection(avatar="",icon("line_weight", color = "blue")),
                            itemsection(slider(0.1:0.1:10, :line_width, label =  false, color = "blue",disable=:enable2DControls))
                            ])
                            ,
                            item([
                            # p("Line Width:"),
                            itemsection(avatar="",icon("format_list_numbered_rtl", color = "blue")),
                            
                            itemsection(slider(2:0.1:25, :num_contours, label = false, color = "blue",disable=:enable2DControls))
                            ])
                            ,
                            item([
                            # p("Line Width:"),
                            itemsection(avatar="",icon("gradient", color = "blue")),
                            itemsection(slider(0:0.01:0.99, :mid_point, label = false, color = "blue",disable=:enable2DControls))
                            ]),
                            item([
                            toggle(
                                :toggleCentreLabel,
                                :centreBoolean,
                                leftlabel = true,
                                # label= :toggleLabel,
                                color = "blue",
                                var"true-value" = "Yes",
                                var"false-value" = "No",
                                disable=:enable2DControls)
                            ]),
                            item([
                                range(0:5:100, :Range_y, reverse=true, vertical = true, labelalways = true,markers = true)
                            ]),
                            item([range(0:5:100, :Range_x, label = true, labelalways = true,markers = true)]),   
                            item([
                                p("Please note these axes numbers are scaled strangely. There are 100 points plotted between a z-score of -5 and +5. Therefore the mean (0) is at 50. 1 standard deviation below the mean (usually labelled -1) is at 40. 1 standard deviation above the mean (usually labelled as +1) is at 60. You should update these axes labels in illustrator when you have found a plot you like.")])
    
                        ])

                    ]),
                ],
            ),
        ]),
            cell(class="col-md-6", [
                plot(:trace, layout=:thisLayout),
                # Plots.plot(:trace),
                plot(:trace2,layout=:thisLayout2)
                ])
        ])
 
    
    
end

@page("/", ui)
end


