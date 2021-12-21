push!(LOAD_PATH, MDPpath * "Preprocessing");
push!(LOAD_PATH, MDPpath * "Potential");

# Set Julia's PATH enviroment variable so that MDP can call external programs
ENV["PATH"] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";

using Preprocessing, Potential

# load Snap library 
Potential.loadsnap(MDPpath)






