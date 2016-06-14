#Pkg.add("DataFrames")
#Pkg.add("MAT")
using DataFrames
using MAT

using HDF5;

function gc_sac_contacts(blocksize = [15; 11])

blocksize = vec(blocksize)

h5path = "/usr/people/smu/seungmount/research/Alex/Skeleton/e2198_warped.h5";
sacsomaCSV = "/data/smu/dev/e2198_Ca_imaging/data/sac_soma_centers_m2_warped.csv"
soma_low_cut_off = 651;
soma_high_cut_off = 1012;

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

#=
somaDict = [
70244 518 609 1366
70243 523 4725 2237
]
somaDict = somaDict[:, [1; xyz[1:2]+1]];
somaDict = Dict((somaDict[_,1], somaDict[_,2:end]) for _ in 1:size(somaDict,1));
=#

df = readtable(sacsomaCSV, header = false)
somaDict = Dict((df[_,1], vec(Array(df[_,2:3]))) for _ in 1:size(df,1));
#display(somaDict);
#return
#somaDict = map!((x) -> (x[1], x[2][xyz[1:2]]), somaDict);


h5fid = h5open(h5path, "r");
voldata = h5fid["/main"];
#% dims = size(data)  #(2048,5376,3456)
voldim = [size(voldata)...];
griddim = ceil(Int, (voldim[xyz[1:2]]-2) ./ blocksize);

#grid_GCs = Array{Set}(griddim...);
#T = eltype(voldata);
grid_GCs = [Set{eltype(voldata)}() for ii in 1:griddim[1], jj in 1:griddim[2]];

grid_denom = zeros(Int, 360, griddim...);
display(size(grid_denom))

@time if true
gc_num = Dict(0 => zeros(Int, 360));

#for gridy = 1:100 #:griddim[2]
#	for gridx = 1:griddim[1]
for gridy = 1:griddim[2]
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
		#gclist = Set{typeof(blockdata[1])}();
		#grid_GCs[gridx, gridy] = gclist;
		gclist = grid_GCs[gridx, gridy];

		#for z = soma_low_cut_off+1:soma_high_cut_off-1	# has potential of cutting off GC voxels but we probably don't care
		for z = 2 : soma_high_cut_off-soma_low_cut_off	# has potential of cutting off GC voxels but we probably don't care
			for y = 2:jmax-1
				for x = 2:imax-1
					cur = blockdata[x, y, z];
					if cur == 0
						continue
					end
					if cur != blockdata[x-1, y, z] || cur != blockdata[x, y-1, z] || cur != blockdata[x, y, z-1] ||
						cur != blockdata[x+1, y, z] || cur != blockdata[x, y+1, z] || cur != blockdata[x, y, z+1]
						# is surface point
						soma = get(somaDict, cur, []);
						if isempty(soma)   # non-SAC
							push!(gclist, cur);
						else  # SAC
							soma -= grid_offset; 	# might be faster to precompute these
							vec = ([x; y] - soma) .* resolution_xy;
							angle = round(Int, rad2deg( atan2(vec[2], vec[1]) ) );
							angle = angle<1 ? angle+360 : angle
							#grid_denom(angle, gridx, gridy) += 1; # WTH: this creates a generic function....
							grid_denom[angle, gridx, gridy] += 1;

							# numerator
							contactinggcs = unique([
								 blockdata[x-1, y, z],  blockdata[x, y-1, z], blockdata[x, y, z-1],
								 blockdata[x+1, y, z],  blockdata[x, y+1, z], blockdata[x, y, z+1]]);
							for gc in contactinggcs
								get!(gc_num, gc, zeros(Int, 360))[angle] += 1;
							end
						end
					end
				end
			end
#end	# @time
		end
	end
end
# what an ugly tail of nested loops


end # @time
@time if true
display(length(grid_GCs))
#allGCs = union(grid_GCs...);	#oom kill
allGCs = grid_GCs[1];
for gclist in grid_GCs
	union!(allGCs, gclist);
end
print("n GCs = ", length(allGCs));

end # @time
@time if true
gc_denom = Dict( (_, zeros(Int, 360)) for _ in allGCs );
end # @time
@time if true
for (gclist, denom) in zip(grid_GCs, grid_denom[:, _] for _ in 1:length(grid_GCs))
	for gc in gclist
		gc_denom[gc] += denom;
	end
end
end # @time

#=
# doens't work
# convert dict to array to save as matlab cell array 
gc_denom_mat = Array{Nullable{valtype(gc_denom)}}(99999);
gc_denom_mat = Array{valtype(gc_denom)}(99999);
gc_denom_mat = [Nullable{valtype(gc_denom)}() for _ = 1:99999];
gc_denom_mat = Array{Any}(99999);

gc_num_mat = Array{valtype(gc_num)}(99999);
for (id,val) in gc_denom
	gc_denom_mat[id] = val;
end
for (id,val) in gc_num
	gc_num_mat[id] = val;
end
=#
gc_denom_mat = hcat([[elem...] for elem in gc_denom]...);
gc_num_mat = hcat([[elem...] for elem in gc_num]...);

#file = matopen("gc_sac_contacts.mat", "w")
matwrite("gc_sac_contacts.mat", Dict("gc_denom" => gc_denom_mat, "gc_num" => gc_num_mat));

return gc_num, gc_denom

end # function

#=
for chunk / grid point
	for voxel
		if not in on/off layer (maybe not a problem - soma already cut off by Jinseop [how about SAC soma?])
			skip
		end
		if gc [actually omni has GC list already, but probably not in jinseop's data]
			grid_GCs(current location) += gc_id;
		elseif sac_surface
			angle = angleToSACsoma;
			grid_denom(current location)(angle) += 1;	% should count in the resolution...? maybe not... actual contacts are likely not
			OR:
			keep a point list for points for each SAC, and summarize at end of chunk.
			add_point_to_list_of_SAC(sac_id), and compute SAC angle 
		end
	end
end
initialize gc_denom(gc_id);
for grid point
	for gcs in grid
		gc_denom(gc) += grid_denom(current location);
	end
end
=#

#=
verified: 17161 vs 70244

for gridy = 1:100
	for gridx = 1:griddim[1]

julia> find(rr[2][17161])
4-element Array{Int64,1}:
 318
 319
 329
 334

julia> rr[2][17161][find(rr[2][17161])]
4-element Array{Int64,1}:
  2
  2
 15
  1

julia> rr[1][17161][find(rr[2][17161])]
4-element Array{Int64,1}:
  556
  613
 1580
 1759

=#
