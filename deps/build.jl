# deps/build.jl
const BIN_DIR = joinpath(homedir(), ".julia", "bin")
mkpath(BIN_DIR)

# The name of the executable command you want the user to type
const EXE_NAME = "polyorigin" 
const EXE_PATH = joinpath(BIN_DIR, EXE_NAME)

# Path to the actual script inside the installed package
const SCRIPT_PATH = joinpath(@__DIR__, "..", "cli", "cli_polyorigin.jl")

@info "Installing CLI shortcut to $EXE_PATH"

if Sys.iswindows()
    # Create a batch file for Windows
    open(EXE_PATH * ".bat", "w") do f
        write(f, "@echo off\njulia --project=@. \"$SCRIPT_PATH\" %*")
    end
else
    # Create a bash wrapper script for macOS/Linux
    open(EXE_PATH, "w") do f
        write(f, "#!/bin/sh\nexec julia --project=@. \"$SCRIPT_PATH\" \"\$@\"")
    end
    chmod(EXE_PATH, 0o755) # Make it executable
end
