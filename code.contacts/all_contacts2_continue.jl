#Pkg.add("DataFrames")
#Pkg.add("MAT")
#Pkg.add("JLD")
using DataFrames
using MAT
using JLD

using HDF5;

function all_contacts2_continue(ndim = 3)


#% light axis was x, reorder light axis to z.
xyz = [2; 3; 1];
ijk = [3; 1; 2];
resolution = [16.5; 16.5; 23];
resolution_xy = resolution[xyz[1:2]];

h5path = "/usr/people/smu/seungmount/research/Alex/retina/e2198_warped.h5";
h5fid = h5open(h5path, "r");
voldata = h5fid["/main"];
#% dims = size(data)  #(2048,5376,3456)
voldim = [size(voldata)...];


exceptions = [];


@time begin
display("reading raw contact array")
basepath = "/usr/people/smu/dev/e2198_Ca_imaging/data/contacts/raw3d";
fid = open("$basepath/allcontacts2.jls","r");
contactmap = deserialize(fid);
close(fid);
#deserialize(open("deps/datasets/SNEMI3D/ds_test/affinities.jls","r"))
end # @time



##=
maxl = 0;
maxn = 0;
for y = 2:voldim[xyz[2]]-1
    for x = 2:voldim[xyz[1]]-1
        n = length(contactmap[x,y])
        if n > maxn
            maxn = n;
        end
    end
end
# =#

display(maxl)
display(maxn)   #244

# TODO: maybe put this into the initial pass, to get the cells with no contact too?
@time begin
ids = Set{Int32}();
for y = 2:voldim[xyz[2]]-1
    for x = 2:voldim[xyz[1]]-1
        #union!(ids, unique(contactmap[x,y])); #correct for 2d only
        union!(ids, unique(contactmap[x,y][1:ndim:length(contactmap[x,y])]));
    end
end
end # @time
display(length(ids))  #1796




##=
@time begin
display("writing remaining cells")
for id in ids
    if ispath("$basepath/$id.mat")
        #display("skipping $id")
        continue;
    end
    display(id)

    @time begin
    #cellContacts = zeros(Int32, 3, 0);  # cell2, x, y  #2d   # cell2, x, y, z #3d
    cellContacts = zeros(Int32, 0);
    for y = 2:voldim[xyz[2]]-1
        for x = 2:voldim[xyz[1]]-1
            pix = contactmap[x,y];
            for ii = 1:ndim:length(pix)
                if pix[ii] == id
                    #cellContacts = [cellContacts [pix[ii+1]; x; y]];  % too slow, also turns out this promotes the type to int64 rather than int32
                    #append!(cellContacts, [pix[ii+1]; x; y]);
                    append!(cellContacts, [pix[ii+1]; x; y; pix[ii+2]]);    #3d
                end
            end
        end
    end
    #cellContacts = reshape(cellContacts, (3, div(length(cellContacts),3))); # WTH: must use tuple for dims
    cellContacts = reshape(cellContacts, (1+ndim, div(length(cellContacts),4))); # WTH: must use tuple for dims  #3d
    end # @time
    save("$basepath/$id.jld", Dict("contacts" => cellContacts),
        compress=true);
    matwrite("$basepath/$id.mat", Dict("contacts" => cellContacts));
end
end # @time
# =#

return  #3d

#2d
contactmap2 = zeros(Int32, voldim[xyz[1:2]]..., maxn);
    for y = 2:voldim[xyz[2]]-1
        for x = 2:voldim[xyz[1]]-1
            pix = contactmap[x,y];
            contactmap2[x, y, 1:length(pix)] = pix;
        end
    end

#return contactmap, exceptions

# Use serialize() to save as JLS?
# Each element in the array is saved as its own HDF5 dataset, infeasible due to overhead.
# Also see https://github.com/JuliaLang/julia/issues/7909
# eg. 
# > whos()
#  b     78 KB     100×100 Array{Array{Int32,1},2}
# > save("all_contacts3.jld", Dict("x" => b,), compress=true);
# A 2.97MB file is saved. - roughly 300B overhead per element.
# And this is not done saving after 5GB on disk:   confusion_mat2 145654 KB     5376×3456 Array{Array{Int32,1},2}


matwrite("$basepath/all_contacts2.mat", Dict(
    "contactmap" => contactmap2,
    ));
save("$basepath/all_contacts2.jld", Dict("contactmap" => contactmap2,
    "exceptions" => exceptions,
    ), compress=true);

return contactmap2, exceptions

end # function

# ==#
