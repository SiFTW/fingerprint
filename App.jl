module App
using Main.StatisticAnalysis
using GenieFramework
using PlotlyBase
# using PlotlyJS
# using Plots
using Statistics
using KernelDensity
using Interpolations
using DataFrames
using Colors
using CSV
# Plots.plotly()
@genietools



@app begin
    
    @out fileIndex=1
    @out dataType="stained"
    @in xlabel =""
    @in ylabel =""
    @out p="upload first sample/cell line"
    const FILE_PATH = joinpath("public", "uploads")
    mkpath(FILE_PATH)
    @out upfiles = readdir(FILE_PATH)
    @out thisLayout = PlotlyBase.Layout(title="Fingerprint")
    @out thisLayout2 = PlotlyBase.Layout(title="",yaxis_visible=false,xaxis_visible=false,showlegend = true,height=200,plot_bgcolor = "white")

    @private data = DataFrame()
    @in selected_file = ""
    @in selected_column1 = ""
    @in selected_column2 = ""
    @out columns = []
    @out uploadBoolean=true
    @out highlightCol1=false
    @out highlightCol2=false
    @out highlightLabel1=true
    @out highlightLabel2=false
    @out hideSelections1=false
    @out hideSelections2=false
    @out columnSelectText1 = ""
    @out columnSelectText2 = ""
    @out instructionText = ""
    @out columnDict = Dict()
    @out lastStained=""
    @out stainedUnstainedDict = Dict()
    @out allDatasets = []
    @out namingBoolean = false
    @out columnBoolean = true
    @out fileUploadLabel = ""
    @in dimensions = "2D"
    @out toggleLabel = "Dimensions: 2D"
    @out normalizeLabel = "Normalise: No"
    @in normalizeHeight = "No"
    @in line_width= 4
    @in num_contours=5
    @in mid_point=0.1
    @out colorSchemesDict = Dict("colorblind" => :seaborn_colorblind, "bright" => :seaborn_bright)
    @out colorSchemesKeys = ["colorblind","bright"]
    # p=Plots.plot([1,2,3],show=false)
    # @out trace = [Plots.scatter(x=[],y=[])]
    @out trace = []
    @out trace2 = []
    # @out trace = []
    @in enable2DControls=false
    @in colors=["#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","gray"]
    @in Range_x = RangeData(30:60)
    @in Range_y = RangeData(30:60)
    @in customHex = "#CC79A7"
    @in lineToEdit = "select experiment"
    @in colorBoolean = true
    @in linecolorDict=Dict()
    @in delete_line = false
    @out toggleCentreLabel = "Add centre dot: No"
    @in centreBoolean = "No"

    function updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)

        # collect the stained and unstained file names for 
        # mergedDF=DataFrame(repeat([],1,3),["experiment",xlabel,ylabel])
        mergedDF=DataFrame([String[],Float64[],Float64[]], ["experiment", xlabel,ylabel])
        allStained=collect(keys(stainedUnstainedDict))
        for i in 1:length(allStained)
    
            thisStained=allStained[i]
            thisUnstained=stainedUnstainedDict[allStained[i]]
            dataStained = CSV.read(joinpath(FILE_PATH, thisStained), DataFrame)
            dataUnstained = CSV.read(joinpath(FILE_PATH, thisUnstained), DataFrame)
    
            numberOfEvents=size(dataStained,1)
            # ArrayToPopulate=zeros(numberOfEvents,3)
            ArrayToPopulate=Array{Any}(undef, numberOfEvents, 3)
            (column1Stained,column2Stained)=columnDict[thisStained]
            (column1UnStained,column2UnStained)=columnDict[thisStained]
    
            eachMeasurement1=dataStained[:,column1Stained]
            eachMeasurement2=dataStained[:,column2Stained]
            eachMeasurementUnStained1=dataStained[:,column1UnStained]
            eachMeasurementUnstained2=dataStained[:,column2UnStained]
            thisUnstainedMFI1=median(eachMeasurementUnStained1)
            thisUnstainedMFI2=median(eachMeasurementUnstained2)
            ToAdd1=eachMeasurement1.-thisUnstainedMFI1
            ToAdd2=eachMeasurement2.-thisUnstainedMFI2
            ArrayToPopulate[:,1].=thisStained
            ArrayToPopulate[:,2].=Float64.(ToAdd1)
            ArrayToPopulate[:,3].=Float64.(ToAdd2)
            ArrayToPopulate=DataFrame(ArrayToPopulate,names(mergedDF))
            ArrayToPopulate=filter(row -> (sum([row[2],row[3]])>2),  ArrayToPopulate)
            append!(mergedDF,ArrayToPopulate)
        end
        
        CSV.write("testOutput.csv", mergedDF)
        mergedDF[:,2].= (mergedDF[:,2] .- mean(mergedDF[:,2])) ./ std(mergedDF[:,2])
        mergedDF[:,3].= (mergedDF[:,3] .- mean(mergedDF[:,3])) ./ std(mergedDF[:,3])
        CSV.write("testOutputStd.csv", mergedDF)
        trace=[]
        trace2=[]
        thisLayout=PlotlyBase.Layout()
        thisLayout2=PlotlyBase.Layout()
        for i in 1:length(allStained)
            thisDF=mergedDF[(mergedDF.experiment .== allStained[i]), :]
            k = kde((thisDF[:,2], thisDF[:,3]),npoints=(100,100),boundary=((-5,5),(-5,5)))
            if normalizeHeight == "Yes"
              k.density=k.density./maximum(k.density)
              
            end
            # 

            if dimensions == "3D"
                # k.density=k.density./num_contours
                
                thisTrace = surface(z=k.density, line_width=line_width,opacity=0.3,showscale=false,colorscale=[[0, "white"], [1, linecolorDict[allStained[i]]]],contours_z=attr(
                    show=true,
                    usecolormap=false,
                    color=linecolorDict[allStained[i]],
                    highlightcolor="limegreen",
                    project_z=true,
                    line_width=line_width,
                    z=attr(show=true, start= 0.1, size= 0.1),
                    z_end=0.5
                ))
                # plot!(p1,k,levels=[0.1,0.2,0.3,0.4,0.5],label=allStained[i],ylims=(-3,3),xlim=(-3,3),linewidth = 3)
                # push!(trace,thisTrace)
                trace=vcat(trace,[thisTrace])
                thisLayout = PlotlyBase.Layout(title=xlabel*" against "*ylabel*" fingerprint",
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
                thisTrace = contour(;z=k.density,line_width=line_width,colorscale=[[0, "white"],[mid_point,linecolorDict[allStained[i]]],[1, linecolorDict[allStained[i]]]], contours_coloring="lines",showscale=false,ncontours=num_contours)
                
                # contours_z=attr(z=attr(show=true, start= 0.1, size= 0.1),
                # z_end=0.5))
                
                trace=vcat(trace,[thisTrace])
                if centreBoolean =="Yes"
                    yCentre=median(thisDF[:,2])
                    xCentre=median(thisDF[:,3])
                    xCentre=50+(10*xCentre)
                    yCentre=50+(10*yCentre)
                    extraPoint=scatter(x=[xCentre,xCentre],y=[yCentre,yCentre],mode="markers",color=linecolorDict[allStained[i]],marker=attr(size=12,color=linecolorDict[allStained[i]]))
                    trace=vcat(trace,[extraPoint])
                end
                thisLayout = PlotlyBase.Layout(title=xlabel*" against "*ylabel*" fingerprint",
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
            thisTrace2=scatter(x=[1,2],y=[1,1],size=(2,2),mode="lines",line_width=line_width,line=attr(color=linecolorDict[allStained[i]],width=line_width),name=allStained[i])
            trace2=vcat(trace2,[thisTrace2])
            thisLayout2=PlotlyBase.Layout(title="legend",yaxis_visible=false,xaxis_visible=false,showlegend = true,height=200,plot_bgcolor = "white")
                
    
        end
        return trace, thisLayout,trace2,thisLayout2
    end

    @onchange dimensions begin    
        toggleLabel= "Dimensions: " *dimensions
        trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)
        enable2DControls=(dimensions=="3D")
    end

    @onchange centreBoolean begin
        if centreBoolean == "Yes"
            toggleCentreLabel= "Add centre dot: Yes"
        else
            toggleCentreLabel= "Add centre dot: No"
        end
        trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)
    end

    @onchange lineToEdit begin    
        colorBoolean=false
    end

    @onchange customHex begin
        linecolorDict[lineToEdit]=customHex
        trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)

    end

    @onbutton delete_line begin
        delete!(stainedUnstainedDict,lineToEdit)
        trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)

    end

    @onchange line_width, num_contours, mid_point, Range_x, Range_y begin    
        trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)
    end
    
    @onchange normalizeHeight begin
        normalizeLabel= " Normalise: " * normalizeHeight
        trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)
    end

    @in tab_selected = "Upload"
    @onchange selected_column1 begin
        highlightCol1 = false
        highlightCol2 = true
        hideSelections1=false
        
    end
    
    
    



    @onchange selected_column2 begin
        hideSelections2=false
        columnDict[selected_file]=(selected_column1,selected_column2)
        columnBoolean=true
        highlightCol2 = false
        if dataType == "stained"
            dataType="unstained"
            notify(__model__,"Now upload control for this sample")
            instructionText = "Return to file upload and upload the control experiment matching this sample."
            uploadBoolean = false
            lastStained=selected_file
            fileUploadLabel="Upload control CSV file for this dataset:"

        else
            dataType="stained"
            notify(__model__,"Now upload next sample if there is one")
            fileIndex=fileIndex+1
            instructionText = "Return to file upload and upload the next sample if there is one."
            fileUploadLabel="Upload fully stained CSV file for the next dataset:"

            uploadBoolean = false
            # save which file is the stained file and which is the unstained pairing in the dictionary
            stainedUnstainedDict[lastStained]=selected_file
            allDatasets=collect(keys(stainedUnstainedDict))
            
            # at this point we have all the data we need to make or update a fingerprint
            
            linecolorDict[lastStained]=colors[length(linecolorDict)+1]

            trace, thisLayout,trace2,thisLayout2 = updatePlot(xlabel,ylabel,stainedUnstainedDict,FILE_PATH,columnDict,normalizeHeight,dimensions,line_width,num_contours,mid_point,linecolorDict,Range_x,Range_y,centreBoolean)

            
        end
    end
    @onchange selected_file begin
        data = CSV.read(joinpath(FILE_PATH, selected_file), DataFrame)
        columns = names(data)
        uploadBoolean=true
        namingBoolean = true
        columnBoolean = false
        hideSelections1=true
        hideSelections2=true
    end

    @onchange xlabel,ylabel,fileuploads begin
        # trace = [scatter()]
        layout = PlotlyBase.Layout(title=xlabel*" against "*ylabel*"",yaxis_title_text=ylabel, xaxis_title_text=xlabel,height=800)
    end

    @onchange xlabel begin
        columnSelectText1="select column for:"*xlabel
        highlightLabel1=false
        highlightLabel2=true
    end
    @onchange ylabel begin
        columnSelectText2="select column for:"*ylabel
        uploadBoolean = false
        fileUploadLabel = "Upload fully stained CSV file for this dataset:"
        highlightLabel1=false
        highlightLabel2=false
    end

    @onchange fileuploads begin
        if ! isempty(fileuploads)
            @info "File was uploaded: " fileuploads
            filename = fileuploads["name"]
            

            try
                isdir(FILE_PATH) || mkpath(FILE_PATH)
                mv(fileuploads["path"], joinpath(FILE_PATH, string(fileIndex)*"_"*dataType*"_"*filename), force=true)
            catch e
                @error "Error processing file: $e"
                notify(__model__,"Error processing file: $(fileuploads["name"])")
            end
            selected_file=string(fileIndex)*"_"*dataType*"_"*filename
            fileuploads = Dict{AbstractString,AbstractString}()
        end
        upfiles = readdir(FILE_PATH)
    end
    @event uploaded begin
        @info "uploaded"
        notify(__model__, "File was uploaded. Now select column.")
        columnBoolean = false
        highlightCol1 = true
        fileUploadLabel=""
        uploadBoolean=true
    end
    @event rejected begin
        @info "rejected"
        notify(__model__, "Please upload a valid file")
    end
    # @in N = []
    # @in M = []
    # @out Z = 0.0
    # @out trace = [scatter(x=[],y=[])]
    # @out layout = PlotlyBase.Layout(title="Histogram plot")

    # @onchange M begin
    #     xpoints = split(M[1],",")
    #     ypoints = split(N[1],",")
    #     k = kde((xpoints, ypoints),npoints=(res,res),boundary=((-5,10),(-5,10)))
    #     k.density=k.density./maximum(k.density)
    #     z_data=k
    #     trace = [scatter()]
    #     Z=length(xpoints)
    # end
    
end



function ui()
    # plotlyjs()
    # Plots.plotly()
   
    row([
        cell(class="col-md-5", [
            
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
                        p("Upload your data, one experiment at a time."),
                        textfield("First protein name:", :xlabel,placeholder="Bcl2",disable=:namingBoolean,standout=:highlightLabel1),
                        textfield("Second protein name:", :ylabel,placeholder="Bclxl",disable=:namingBoolean,standout=:highlightLabel2),
                        # p("upload files"),
                        # p("{{p}}", class="st-module"),
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
                        list(
                            bordered = true,
                            separator = true,
                            [
                                p("Select a line to edit"),
                                item(
                                    Stipple.select(:lineToEdit; options=:allDatasets,disable=false,usechips=true,clearable=true)
                                ),
                                item(
                                    textfield("Color (as hex code):", :customHex,placeholder="#CC79A7",disable=:colorBoolean)
                                ), 
                                item(
                                    btn("Delete", icon = "delete", color = "red", class = "q-mr-sm",@click(:delete_line)),
                                ),                                
                            ]
                        )
                    ]),

                
                    tabpanel(name = "Tune", [
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
                            itemsection(slider(0.1:0.1:10, :line_width, label = "Line Width", color = "blue",disable=:enable2DControls))
                            ])
                            ,
                            item([
                            # p("Line Width:"),
                            itemsection(avatar="",icon("format_list_numbered_rtl", color = "blue")),
                            
                            itemsection(slider(2:0.1:25, :num_contours, label = "num_contours", color = "blue",disable=:enable2DControls))
                            ])
                            ,
                            item([
                            # p("Line Width:"),
                            itemsection(avatar="",icon("gradient", color = "blue")),
                            itemsection(slider(0:0.01:0.99, :mid_point, label = "mid_point", color = "blue",disable=:enable2DControls))
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
 
    # row([
    #      cell(class="col-md-3", [
    #         # p("Plot Settings"),
    #         textfield("Protein #1", :xlabel,placeholder="Bcl2",disable=:namingBoolean),
    #         textfield("Protin #2", :ylabel,placeholder="Bclxl",disable=:namingBoolean),
    #         # p("upload files"),
    #         # p("{{p}}", class="st-module"),
    #         uploader( multiple = false,
    #             accept = ".csv",
    #             maxfilesize = 1024*1024*150, # bytes
    #             maxfiles = 30,
    #             autoupload = true,
    #             hideuploadbtn = true,
    #             label = "Upload datasets",
    #             nothumbnails = true,
    #             style="max-width: 95%; width: 95%; margin: 0 auto;",
    #             disable=:uploadBoolean,
    #             @on("rejected", :rejected),
    #             @on("uploaded", :uploaded)
    #             ),
    #         p("{{columnSelectText1}}"),
    #         Stipple.select(:selected_column1; options=:columns,disable=:columnBoolean),
            
    #         p("{{columnSelectText2}}"),
    #         Stipple.select(:selected_column2; options=:columns,disable=:columnBoolean),
    #         p("{{instructionText}}")
    #         ]),
    #         cell(class="col-md-8", [
    #             plot(:trace, layout=:layout),
    #             # Plots.plot(:trace),
    #             toggle(
    #                 :toggleLabel,
    #                 :dimensions,
    #                 color = "blue",
    #                 val= 2,
    #                 var"true-value" = "3D",
    #                 var"false-value" = "2D",
    #             )

    #             ])
    #     ])
    
end

@page("/", ui)
end


# module App
# using GenieFramework
# @genietools

# @app begin
#     @out firstX = [1, 2, 3, 4, 5]
#     @out firstY = [5, 0, 3, -1, 2]
#     @out secondX = [1, 2, 3, 4, 5]
#     @out secondY = [1, 4, 6, 4, 4]
# end

# function ui()

#     row([
#         cell(class="st-col col-3", [
#             h1("A simple dashboard"),
#             slider(1:1000, :N),
#             p("The average of {{N}} random numbers is {{m}}", class="st-module"),
#             plot(:trace, layout=:layout)
#         ])
#     ])
# end
    
# @page("/", "plot-in-HTML.html")
# Server.isrunning() || Server.up()
# end