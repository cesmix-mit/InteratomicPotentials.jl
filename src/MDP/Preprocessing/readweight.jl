#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function readweight(app, config)

if app.weightmode == 0
    tmp = ones(1,2*config.nconfigs);
    tmp[1:config.nconfigs] .= app.we;
    tmp[(config.nconfigs+1):2*config.nconfigs] .= app.wf;
elseif app.weightmode == 1 # binary    
    filename = app.weightfile;
    tmp = reinterpret(Float64,read(filename));    
elseif app.weightmode == 2 # text    
    filename = app.weightfile;
    tmp = Main.DelimitedFiles.readdlm(filename);
end

# determine weights
if (app.dftdata == 1)
    config.we = reshape(tmp[1:config.nconfigs],1,config.nconfigs);
    config.wf = reshape([],0,2);
elseif (app.dftdata == 2)
    config.we = reshape([],0,2);
    config.wf = reshape(tmp[1:config.nconfigs],1,config.nconfigs);
elseif (app.dftdata == 3)
    config.we = reshape(tmp[1:config.nconfigs],1,config.nconfigs);
    config.wf = reshape(tmp[(config.nconfigs+1):2*config.nconfigs],1,config.nconfigs);
end

return config

end
