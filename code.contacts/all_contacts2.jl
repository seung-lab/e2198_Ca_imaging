#Pkg.add("DataFrames")
#Pkg.add("MAT")
#Pkg.add("JLD")
using DataFrames
using MAT
using JLD

using HDF5;

function all_contacts2(ndim=3, blocksize = [15; 11], resume::Bool=false) #... now it takes 17hours...no longer io bound?

blocksize = vec(blocksize)

h5path = "/usr/people/smu/seungmount/research/Alex/retina/e2198_warped.h5";
basepath = "/usr/people/smu/dev/e2198_Ca_imaging/data/contacts/raw";
if ndim==3
basepath = "/usr/people/smu/dev/e2198_Ca_imaging/data/contacts/raw3d";
end

#sacsomaCSV = "/data/smu/dev/e2198_Ca_imaging/data/sac_soma_centers_m2_warped.csv"
soma_low_cut_off = 600; #%90 #651;
soma_low_cut_off = 445; #%120
soma_high_cut_off = 1167; #1012;
#on_off_threshold = 831.5;

#=
soma_low_cut_off = 200;
soma_high_cut_off = 1500;
# =#

#=
x=919 is 28% 
x=744 is 62% 
simple math gives
x=445 is 120.091... %
x=1167 is -20.183... %

10 = 1012
80 = 651
=#

#% light axis was x, reorder light axis to z.
xyz = [2; 3; 1];
ijk = [3; 1; 2];
resolution = [16.5; 16.5; 23];
resolution_xy = resolution[xyz[1:2]];



h5fid = h5open(h5path, "r");
voldata = h5fid["/main"];
#% dims = size(data)  #(2048,5376,3456)
voldim = [size(voldata)...];

griddim = ceil(Int, (voldim[xyz[1:2]]-2) ./ blocksize);

#confusion_mat = Dict{UInt64, UInt}();
#contactmap = zeros(Int32, voldim[xyz[1:2]]..., 100);
#contactmap = Array{Dict{UInt64, Int8}, voldim[xyz[1:2]]...);
#contactmap = [Dict{UInt64, Int8}() for x in 1:voldim[xyz[1]], y in 1:voldim[xyz[2]]];
contactmap = [Array{Int32, 1}() for x in 1:voldim[xyz[1]], y in 1:voldim[xyz[2]]];
display(sizeof(contactmap))


if resume    #post processing

    @time begin
    display("reading raw contact array")
    fid = open("$basepath/allcontacts2.jls","r");
    contactmap = deserialize(fid);
    close(fid);
    end # @time

else

@time begin

exceptions = [];
surfacearea = Dict{Int32, Int}();   # note this is in reality voxel count without account to anisotropy
surface_at_grid = Dict{Int32, Array{Int32,2}}(); # Dict(0 => zeros(Int32, griddim...));
confusion_mat = Dict{UInt64, UInt}();

#for gridy = 34:34 #:griddim[2]
#   for gridx = 34:34 #:griddim[1]
for gridy = 1:griddim[2]
    display(gridy)
    for gridx = 1:griddim[1]
        grid_offset = ([gridx; gridy] - 1) .* blocksize;
        grid_max = [gridx; gridy] .* blocksize + 2;
#display(grid_max)
#display(voldim)
        grid_max = min(grid_max, voldim[xyz[1:2]]);  # note: has to be array (not tuple) to take element-wise min
#@time if true

        datarange = (grid_offset[1]+1:grid_max[1], grid_offset[2]+1:grid_max[2], soma_low_cut_off:soma_high_cut_off);
        datarange = datarange[ijk];

        #voldata = h5read(path, "img", slice)
        blockdata = voldata[datarange...];

        blockdata = permutedims(blockdata, xyz);

        imax, jmax = grid_max - grid_offset;
#end # @time
#@time if true

        #for z = soma_low_cut_off+1:soma_high_cut_off-1 # has potential of cutting off GC voxels but we probably don't care
        for z = 2 : soma_high_cut_off-soma_low_cut_off  # has potential of cutting off GC voxels but we probably don't care
            
            for y = 2:jmax-1
                for x = 2:imax-1
                    cur = blockdata[x, y, z];
                    if cur == 0
                        continue
                    end
                    thisKey = UInt64(cur) << 32;
                    coord = grid_offset + [x; y];
                    #coord = [grid_offset + [x; y]; soma_low_cut_off-1+z]; #3d
                    #d = contactmap[coord...];
                    issurface = false;
                    iscontact = false;  # TODO: count total contact (which could be different from counting the same voxel for multiple contact targets)
                    for neigbor in unique([
                                     blockdata[x-1, y, z],  blockdata[x, y-1, z], blockdata[x, y, z-1],
                                     blockdata[x+1, y, z],  blockdata[x, y+1, z], blockdata[x, y, z+1]]);
        
                        if cur != neigbor
                            issurface = true;
                            key = thisKey + neigbor;
                            if neigbor != 0
                                if false    # skip building contactmap
                                if ndim==2  #2d
                                    append!(contactmap[coord...], [cur, neigbor]);
                                else    #3d
                                    append!(contactmap[coord...], [cur, neigbor, soma_low_cut_off-1+z]); #3d
                                end
                                end

                                confusion_mat[key] = get(confusion_mat, key, 0) + 1;
                            end
                        end
                    end
                    if issurface
                        surfacearea[cur] = get!(surfacearea, cur, 0) + 1;

                        get!(surface_at_grid, cur, zeros(Int32, griddim...))[gridx,gridy] += 1;
                    end
                end
            end
#end    # @time
        end
    end
end
# what an ugly tail of nested loops


end # @time

end # if resume==1

display(sizeof(contactmap))

#=
maxl = 0;
maxn = 0;
for y = 2:voldim[xyz[2]]-1
    for x = 2:voldim[xyz[1]]-1
        d = contactmap[x,y];
        l = length(d);
        if l > maxl
            maxl = l;
        end
        n = sum(values(d));
        if n > maxn
            maxn = n;
        end
        if any([values(d)...].<=0)
            warn("overflown")
        end
    end
end
# =#

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
display(length(ids))  #1796 / 1797



@time if !resume
display("writing raw contact array")
fid = open("$basepath/allcontacts2.jls","w+");
serialize(fid, contactmap);
close(fid);
surfacearea_mat = hcat([[key, count] for (key, count) in surfacearea]...);
confusion_mat = hcat([[key >> 32, key << 32 >> 32, count] for (key, count) in confusion_mat]...);

# this is similar to the older all_contacts.mat
matwrite("$basepath/surface_contingency.mat", Dict(
    "confusion_mat" => confusion_mat,
    ));
save("$basepath/surface_contingency.jld", Dict(
    "confusion_mat" => confusion_mat,
    ));
matwrite("$basepath/surface_grid.mat", Dict(
    "surface_at_grid_keys" => collect(keys(surface_at_grid)),
    "surface_at_grid_vals" => collect(values(surface_at_grid)),
    #"surfacearea" => surfacearea_mat,
    ));
save("$basepath/surface_grid.jld", Dict(
    "surface_at_grid" => surface_at_grid,
    ), compress=true);
matwrite("$basepath/surface_area.mat", Dict(
    "surfacearea" => surfacearea_mat,
    ));
save("$basepath/surface_area.jld", Dict(
    "surfacearea" => surfacearea,
    ));
else
    display(isequal(ids), Set{Int32}(keys(surfacearea)))
    display(length(surfacearea))
end
return

##=
@time begin
display("writing each cell")
if resume
    display("writing remaining cells")
end
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
    cellContacts = reshape(cellContacts, (1+ndim, div(length(cellContacts),1+ndim))); # WTH: must use tuple for dims  #3d
    end # @time
    save("$basepath/$id.jld", Dict("contacts" => cellContacts),
        compress=true);
    matwrite("$basepath/$id.mat", Dict("contacts" => cellContacts));
end
end # @time
# =#

if ndim==3
return surfacearea, exceptions
end

#2d only
contactmap2 = zeros(Int32, voldim[xyz[1:2]]..., maxn);
    for y = 2:voldim[xyz[2]]-1
        for x = 2:voldim[xyz[1]]-1
            contactmap2[x, y, 1:length(pix)] = contactmap[x,y];
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


matwrite("all_contacts2.mat", Dict(
	"contactmap" => contactmap2,
    "surfacearea" => surfacearea,
	));
save("all_contacts2.jld", Dict(
    "contactmap" => contactmap2,
    "surfacearea" => surfacearea,
	"exceptions" => exceptions,
    ), compress=true);

return contactmap2, exceptions

end # function

# ==#