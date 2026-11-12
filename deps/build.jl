# deps/build.jl
depot = isempty(DEPOT_PATH) ? joinpath(homedir(), ".julia") : first(DEPOT_PATH)
BIN_DIR = joinpath(depot, "bin")
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

# test if BIN_DIR is in PATH
path_dirs = split(get(ENV, "PATH", ""), Sys.iswindows() ? ';' : ':')

# Normalize paths (handles trailing slashes and Windows drive letter casing)
norm_bin = normpath(BIN_DIR)
in_path = any(p -> !isempty(p) && normpath(p) == norm_bin, path_dirs)

if !in_path 
    msg = if Sys.iswindows() 
        """
        $(BIN_DIR) is not in your PATH.
        To run $EXE_NAME from anywhere in your terminal, add it to your PATH permanently:
            [Environment]::SetEnvironmentVariable('Path', [Environment]::GetEnvironmentVariable('Path', 'User') + ';$(BIN_DIR)', 'User')            
        """
    else
        """
        $(BIN_DIR) is not in your PATH.
        To run 'mypkg' from anywhere, add this line to your ~/.bashrc or ~/.zshrc:
          export PATH="$(BIN_DIR):\$PATH"
        """
    end

    # 1. Write to build.log
    @warn msg

    # 2. Force print directly to the user's terminal
    tty = Sys.iswindows() ? "CON" : "/dev/tty"
    try
        open(tty, "w") do io
            printstyled(io, "Warning: ", color=:yellow, bold=true)
            println(io, msg)
        end
    catch        
        # Ignore in headless environments (CI/Docker)
    end
end
