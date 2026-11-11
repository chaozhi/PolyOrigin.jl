"""
    get_datafile(filename::String;dataname="polyorigin2_datasets")

gets the file path of `filename`` from `polyorigin2_datasets` artifact.
    
"""
function get_datafile(filename::String; dataname="polyorigin2_datasets")
    # Locate (or automatically download) the unpacked artifact root folder
    root_dir = @artifact_str(dataname)    
    datafile = joinpath(root_dir, filename)    
    if !isfile(datafile)
        @error string("filename must be one of ",get_datafiles())
    end
    datafile
end


"""
    get_datafiles(;dataname="polyorigin2_datasets")

lists datafiles from `polyorigin2_datasets` artifact.
    
"""
function get_datafiles(dataname="polyorigin2_datasets"; join=false)    
    root_dir = @artifact_str(dataname)    
    readdir(root_dir; join)    
end


"""
    get_datadir(;dataname="polyorigin2_datasets")

gets example data directory from the `polyorigin2_datasets` artifact.
    
"""
function get_datadir(dataname="polyorigin2_datasets")    
    @artifact_str(dataname)    
end
