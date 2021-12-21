#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

using Libdl

global snaplib

function loadsnap(mdppath::String)    
    libexist = 0
    if Sys.isapple()
        libpath = mdppath * "Potential/cpuSnap.dylib";
    elseif Sys.isunix()
        libpath = mdppath * "Potential/cpuSnap.so";
    elseif Sys.windows()
        libpath = mdppath * "Potential/cpuSnap.so";
    end
    if isfile(libpath)
        libexist = 1;
    end        
    if libexist==0            
        cdir = pwd();
        cd(mdppath * "Potential")        
        cstr = `g++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuSnap.cpp -o cpusnap.o`                
        run(cstr)
        if Sys.isapple()            
            cstr = `g++ -shared cpusnap.o -o cpuSnap.dylib`            
        else
            cstr = `g++ -shared cpusnap.o -o cpuSnap.so`            
        end 
        run(cstr)
        cd(cdir)                   
    end
    global snaplib = Libdl.dlopen(libpath)
end

function closesnap()
    Libdl.dlclose(snaplib)
end

function bnum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxbcount), Cint, (Cint,), twojmax)
    return count
end

function unum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxucount), Cint, (Cint,), twojmax)
    return count
end

function znum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxzcount), Cint, (Cint,), twojmax)
    return count
end

function cgnum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxcgcount), Cint, (Cint,), twojmax)
    return count
end

function factable()

    facto = [
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600,
        6227020800,
        87178291200,
        1307674368000,
        20922789888000,
        355687428096000,
        6.402373705728e+15,
        1.21645100408832e+17,
        2.43290200817664e+18,
        5.10909421717094e+19,
        1.12400072777761e+21,
        2.5852016738885e+22,
        6.20448401733239e+23,
        1.5511210043331e+25,
        4.03291461126606e+26,
        1.08888694504184e+28,
        3.04888344611714e+29,
        8.8417619937397e+30,
        2.65252859812191e+32,
        8.22283865417792e+33,
        2.63130836933694e+35,
        8.68331761881189e+36,
        2.95232799039604e+38,
        1.03331479663861e+40,
        3.71993326789901e+41,
        1.37637530912263e+43,
        5.23022617466601e+44,
        2.03978820811974e+46,
        8.15915283247898e+47,
        3.34525266131638e+49,
        1.40500611775288e+51,
        6.04152630633738e+52,
        2.65827157478845e+54,
        1.1962222086548e+56,
        5.50262215981209e+57,
        2.58623241511168e+59,
        1.24139155925361e+61,
        6.08281864034268e+62,
        3.04140932017134e+64,
        1.55111875328738e+66,
        8.06581751709439e+67,
        4.27488328406003e+69,
        2.30843697339241e+71,
        1.26964033536583e+73,
        7.10998587804863e+74,
        4.05269195048772e+76,
        2.35056133128288e+78,
        1.3868311854569e+80,
        8.32098711274139e+81,
        5.07580213877225e+83,
        3.14699732603879e+85,
        1.98260831540444e+87,
        1.26886932185884e+89,
        8.24765059208247e+90,
        5.44344939077443e+92,
        3.64711109181887e+94,
        2.48003554243683e+96,
        1.71122452428141e+98,
        1.19785716699699e+100,
        8.50478588567862e+101,
        6.12344583768861e+103,
        4.47011546151268e+105,
        3.30788544151939e+107,
        2.48091408113954e+109,
        1.88549470166605e+111,
        1.45183092028286e+113,
        1.13242811782063e+115,
        8.94618213078297e+116,
        7.15694570462638e+118,
        5.79712602074737e+120,
        4.75364333701284e+122,
        3.94552396972066e+124,
        3.31424013456535e+126,
        2.81710411438055e+128,
        2.42270953836727e+130,
        2.10775729837953e+132,
        1.85482642257398e+134,
        1.65079551609085e+136,
        1.48571596448176e+138,
        1.3520015276784e+140,
        1.24384140546413e+142,
        1.15677250708164e+144,
        1.08736615665674e+146,
        1.03299784882391e+148,
        9.91677934870949e+149,
        9.61927596824821e+151,
        9.42689044888324e+153,
        9.33262154439441e+155,
        9.33262154439441e+157,
        9.42594775983835e+159,
        9.61446671503512e+161,
        9.90290071648618e+163,
        1.02990167451456e+166,
        1.08139675824029e+168,
        1.14628056373471e+170,
        1.22652020319614e+172,
        1.32464181945183e+174,
        1.44385958320249e+176,
        1.58824554152274e+178,
        1.76295255109024e+180,
        1.97450685722107e+182,
        2.23119274865981e+184,
        2.54355973347219e+186,
        2.92509369349301e+188,
        3.3931086844519e+190,
        3.96993716080872e+192,
        4.68452584975429e+194,
        5.5745857612076e+196,
        6.68950291344912e+198,
        8.09429852527344e+200,
        9.8750442008336e+202,
        1.21463043670253e+205,
        1.50614174151114e+207,
        1.88267717688893e+209,
        2.37217324288005e+211,
        3.01266001845766e+213,
        3.8562048236258e+215,
        4.97450422247729e+217,
        6.46685548922047e+219,
        8.47158069087882e+221,
        1.118248651196e+224,
        1.48727070609069e+226,
        1.99294274616152e+228,
        2.69047270731805e+230,
        3.65904288195255e+232,
        5.01288874827499e+234,
        6.91778647261949e+236,
        9.61572319694109e+238,
        1.34620124757175e+241,
        1.89814375907617e+243,
        2.69536413788816e+245,
        3.85437071718007e+247,
        5.5502938327393e+249,
        8.04792605747199e+251,
        1.17499720439091e+254,
        1.72724589045464e+256,
        2.55632391787286e+258,
        3.80892263763057e+260,
        5.71338395644585e+262,
        8.62720977423323e+264,
        1.31133588568345e+267,
        2.00634390509568e+269,
        3.08976961384735e+271,
        4.78914290146339e+273,
        7.47106292628289e+275,
        1.17295687942641e+278,
        1.85327186949373e+280,
        2.94670227249504e+282,
        4.71472363599206e+284,
        7.59070505394721e+286,
        1.22969421873945e+289,
        2.0044015765453e+291,
        3.28721858553429e+293,
        5.42391066613159e+295,
        9.00369170577843e+297,
        1.503616514865e+300,
    ]
    
    facto = Float64.(Array(facto))
    
    return facto
    
end
    
function initsna(param, elemradius, elemweight, snapcoeff=nothing)

    ntypes = Int32(param[1]);  
    twojmax = Int32(param[2]);  
    rcutfac = param[3];
    rfac0 = param[4];
    rmin0 = param[5];
    bzeroflag = Int32(param[6]);
    switchflag = Int32(param[7]);
    quadraticflag = Int32(param[8]);
    chemflag = Int32(param[9]);
    bnormflag = Int32(param[10]);
    wselfallflag = Int32(param[11]);   
    
    sna = initializesna(ntypes, twojmax, chemflag)

    facto = factable()    
    ccall(Libdl.dlsym(snaplib, :cpuInitSna), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint), 
            sna.rootpqarray, sna.cglist, facto, sna.idx_max, sna.idxz, sna.idxz_block, 
            sna.idxb, sna.idxb_block, sna.idxu_block, sna.idxcg_block, sna.twojmax)

    sna.radelem[2:end] = elemradius
    sna.wjelem[2:end] = elemweight    
    sna.map[2:end] = Int32.(Array(1:ntypes))

    if snapcoeff === nothing       
        sna.coeffelem = zeros(Float64, sna.ncoeffall) 
    elseif length(snapcoeff) == sna.ntypes*sna.ncoeffall
        sna.coeffelem = snapcoeff
    else
        error("number of snap coefficients is incompatible with the snap paramaters")
    end

    rcutmax = 0.0
    for ielem = 1:ntypes
        rcutmax = max(2.0*elemradius[ielem]*rcutfac,rcutmax)
    end

    for i = 1:ntypes
        for j = 1:ntypes
            cut = (elemradius[i] + elemradius[j])*rcutfac
            sna.rcutsq[j + (i-1)*(ntypes)] = cut*cut
        end
    end

    wself = 1.0;      
    www = wself*wself*wself;
    for j = 0:twojmax
        if (bnormflag == 1)
            sna.bzero[j+1] = www;
        else
            sna.bzero[j+1] = www*(j+1);
        end
    end

    sna.nperdim = sna.ncoeffall-1;
    sna.bnormflag = bnormflag;
    sna.chemflag = chemflag;    
    sna.quadraticflag = quadraticflag;
    sna.switchflag = switchflag;
    sna.bzeroflag = bzeroflag;
    sna.wselfallflag = wselfallflag;
    sna.wself = wself;
    sna.rmin0 = rmin0;
    sna.rfac0 = rfac0;
    sna.rcutfac = rcutfac;
    sna.rcutmax = rcutmax;        

    return sna
end

function snappotential(x, t, a, b, c, pbc, param, elemradius, elemweight, snapcoeff)

    sna = initsna(param, elemradius, elemweight, snapcoeff)
    eatom, fatom, vatom = snappotential(x, t, a, b, c, pbc, sna)
    return eatom, fatom, vatom 
end

function snappotential(x, t, a, b, c, pbc, sna)

    dim, N = size(x)
    
    #idxcg_max = sna.idxcg_max;
    idxu_max = sna.idxu_max;
    idxb_max = sna.idxb_max;
    idxz_max = sna.idxz_max;    
    twojmax = sna.twojmax;
    #ncoeff = sna.ncoeff;
    ncoeffall = sna.ncoeffall;
    ntypes = sna.ntypes;
    nelem = sna.nelements;    
    # ndoubles = sna.ndoubles;   
    # ntriples = sna.ntriples;   
    # nperdim = sna.nperdim;
    bnormflag = sna.bnormflag;
    chemflag = sna.chemflag;    
    quadraticflag = sna.quadraticflag;
    switchflag = sna.switchflag;    
    bzeroflag = sna.bzeroflag;
    wselfallflag = sna.wselfallflag;        

    map = sna.map;
    idxz = sna.idxz;
    #idxz_block = sna.idxz_block;
    idxb = sna.idxb;
    idxb_block = sna.idxb_block;
    idxu_block = sna.idxu_block;
    idxcg_block = sna.idxcg_block;   
    
    rcutmax = sna.rcutmax
    wself = sna.wself;
    rmin0 = sna.rmin0;
    rfac0 = sna.rfac0;
    rcutfac = sna.rcutfac;
    rcutmax = sna.rcutmax;        
    bzero = sna.bzero;
    rootpqarray = sna.rootpqarray;
    cglist = sna.cglist;
    rcutsq = sna.rcutsq;    
    radelem = sna.radelem;
    wjelem = sna.wjelem; 
    coeffelem = sna.coeffelem;           
    
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:N));   
    atomtype = Int32.(t[:])
    alist = Int32.(alist)
    neighlist = Int32.(neighlist)
    neighnum = Int32.(neighnum)

    #pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcutmax*rcutmax);
    pairlist, pairnum = neighpairlist(y, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)    
    rij, ai, aj, ti, tj = neighpairs(y, pairlist, pairnum, atomtype, ilist, alist);                
    ijnum = pairnum[end]
    
    # offset 1 for C++ code
    ilist = ilist .- Int32(1)
    alist = alist .- Int32(1)
    ai = ai .- Int32(1)
    aj = aj .- Int32(1)

    ulisttot_r = zeros(N*idxu_max*nelem)
    ulisttot_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeUi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        map, ai, ti, tj, twojmax, idxu_max, N, ijnum, switchflag, chemflag)    

    ccall(Libdl.dlsym(snaplib, :cpuAddWself2Ui), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, wself, idxu_block, atomtype, map, wselfallflag, chemflag, 
        idxu_max, nelem, twojmax, N)    

    eatom = zeros(N)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeEi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        eatom, ulisttot_r, ulisttot_i, cglist, bzero, coeffelem, ilist, map, atomtype, 
        idxb, idxcg_block, idxu_block, twojmax, idxb_max, idxu_max, nelem, 
        ncoeffall, bnormflag, bzeroflag, wselfallflag, quadraticflag, N)    
    
    ylist_r = zeros(N*idxu_max*nelem)
    ylist_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeYi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, coeffelem, map, atomtype, idxz, idxb_block, 
        idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelem, ncoeffall, bnormflag, N)    

    fatom = zeros(3*N)       
    vatom = zeros(6*N)            
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeFi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, 
        Cint, Cint, Cint, Cint), 
        fatom, vatom, ylist_r, ylist_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        map, ai, ai, aj, ti, tj, twojmax, idxu_max, N, N, ijnum, switchflag, chemflag)    

    return eatom, fatom, vatom 
end

function snapdescriptors(x, t, a, b, c, pbc, param, elemradius, elemweight)

    sna = initsna(param, elemradius, elemweight)
    bi, bd, bv = snapdescriptors(x, t, a, b, c, pbc, sna)
    return bi, bd, bv 
end

function snapdescriptors(x, t, a, b, c, pbc, sna)

    dim, N = size(x)
    
    #idxcg_max = sna.idxcg_max;
    idxu_max = sna.idxu_max;
    idxb_max = sna.idxb_max;
    idxz_max = sna.idxz_max;    
    twojmax = sna.twojmax;
    ncoeff = sna.ncoeff;
    #ncoeffall = sna.ncoeffall;
    ntypes = sna.ntypes;
    nelem = sna.nelements;    
    ndoubles = sna.ndoubles;   
    ntriples = sna.ntriples;   
    #nperdim = sna.nperdim;
    bnormflag = sna.bnormflag;
    chemflag = sna.chemflag;    
    quadraticflag = sna.quadraticflag;
    switchflag = sna.switchflag;    
    bzeroflag = sna.bzeroflag;
    wselfallflag = sna.wselfallflag;        

    map = sna.map;
    idxz = sna.idxz;
    idxz_block = sna.idxz_block;
    idxb = sna.idxb;
    #idxb_block = sna.idxb_block;
    idxu_block = sna.idxu_block;
    idxcg_block = sna.idxcg_block;   
    
    rcutmax = sna.rcutmax
    wself = sna.wself;
    rmin0 = sna.rmin0;
    rfac0 = sna.rfac0;
    rcutfac = sna.rcutfac;
    rcutmax = sna.rcutmax;        
    bzero = sna.bzero;
    rootpqarray = sna.rootpqarray;
    cglist = sna.cglist;
    rcutsq = sna.rcutsq;    
    radelem = sna.radelem;
    wjelem = sna.wjelem; 
    #coeffelem = sna.coeffelem;                   

    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:N));   
    atomtype = Int32.(t[:])
    alist = Int32.(alist)
    neighlist = Int32.(neighlist)
    neighnum = Int32.(neighnum)
    
    pairlist, pairnum = neighpairlist(y, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)    
    rij, ai, aj, ti, tj = neighpairs(y, pairlist, pairnum, atomtype, ilist, alist);                
    ijnum = pairnum[end]


    # offset 1 for C++ code
    ilist = ilist .- Int32(1)
    alist = alist .- Int32(1)
    ai = ai .- Int32(1)
    aj = aj .- Int32(1)

    n = idxu_max*ijnum
    ulist_r = zeros(n)
    ulist_i = zeros(n)
    dulist_r = zeros(3*n)
    dulist_i = zeros(3*n)    
    ccall(Libdl.dlsym(snaplib, :cpuComputeUij), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, 
        Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
        ulist_r, ulist_i, dulist_r, dulist_i, rootpqarray, rij, wjelem, radelem, rmin0, 
        rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, ijnum, switchflag)    

    ulisttot_r = zeros(N*idxu_max*nelem)
    ulisttot_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuZeroUarraytot2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},  Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, wself, idxu_block, atomtype, map, ai, wselfallflag, chemflag, 
        idxu_max, nelem, twojmax, N)    

    ccall(Libdl.dlsym(snaplib, :cpuAddUarraytot), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, ulist_r, ulist_i, map, ai, tj, idxu_max, N, ijnum, chemflag)    

    zlist_r = zeros(idxz_max*ndoubles*N)
    zlist_i = zeros(idxz_max*ndoubles*N)
    ccall(Libdl.dlsym(snaplib, :cpuComputeZi2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, idxcg_block, twojmax, 
        idxu_max, idxz_max, nelem, bnormflag, N)    

    blist = zeros(idxb_max*ntriples*N)
    ccall(Libdl.dlsym(snaplib, :cpuComputeBi2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, atomtype, 
        map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
        nelem, bzeroflag, wselfallflag, chemflag, N)    
    
    dblist = zeros(idxb_max*ntriples*dim*ijnum)          
    ccall(Libdl.dlsym(snaplib, :cpuComputeDbidrj), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, idxz_block, map, ai, tj, 
        twojmax, idxb_max, idxu_max, idxz_max, nelem, bnormflag, chemflag, N, ijnum)    
    
    bi = zeros(ntypes*ncoeff)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrum), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},  
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint), 
        bi, blist, ilist, atomtype, N, ncoeff, ncoeff, ntypes, quadraticflag)    

    bd = zeros(ntypes*ncoeff*3*N)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrumDeriv), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},  
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        bd, blist, dblist, ai, ai, aj, ti, N, ijnum, ncoeff, ncoeff, ntypes, quadraticflag)    
    
    bv = zeros(6*ntypes*ncoeff)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrumVirial), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},   
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        bv, blist, dblist, rij, ai, ti, N, ijnum, ncoeff, ncoeff, ntypes, quadraticflag)  

    bd = reshape(bd, (3,N,ntypes*ncoeff))
    bv = reshape(bv, (6,ntypes*ncoeff))
    bv = bv[[1; 2; 3; 6; 5; 4],:]

    return bi, bd, bv           
end

