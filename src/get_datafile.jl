"""
    get_datafile(filename::String;dataname="polyorigin2_datasets")

resolves the local cached path for the `polyorigin2_datasets` artifact.
    
"""
function get_datafile(filename::String; dataname="polyorigin2_datasets")
    # Locate (or automatically download) the unpacked artifact root folder
    root_dir = @artifact_str(dataname)    
    datafile = joinpath(root_dir, dataname, filename)    
    datafile
end


"""
    get_datapath(;dataname="polyorigin2_datasets")

resolves the local cached path for the `polyorigin2_datasets` artifact.
    
"""
function get_datadir(dataname="polyorigin2_datasets")    
    root_dir = @artifact_str(dataname)    
    data_dir = joinpath(root_dir, dataname)    
    data_dir
end
