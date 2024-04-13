using DelimitedFiles
path = homedir() * "/Repos/ExternalBJJ/data/fidelity/" 
filenames = readdir(path)
files = path .* filenames
files[2], readdlm(files[2])
filenames

pattern = r"nfid(\d+)([a-zA-Z]+)(\d+)\.dat"
match_result = match(pattern, files[2])
match_result.captures
fileloc(nn, type, u) = path * "nfid" * string(nn) * type * replace(string(u), "." => "") * ".dat" 
fileloc(200, "sta", 0.4)
