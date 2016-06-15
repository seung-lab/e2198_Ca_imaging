#Pkg.add("DataFrames")
#Pkg.add("MAT")
#Pkg.add("JLD")
using DataFrames
using MAT
using JLD

using HDF5;

function gc_sac_contacts(blocksize = [15; 11])

blocksize = vec(blocksize)

h5path = "/usr/people/smu/seungmount/research/Alex/Skeleton/e2198_warped.h5";
sacsomaCSV = "/data/smu/dev/e2198_Ca_imaging/data/sac_soma_centers_m2_warped.csv"
soma_low_cut_off = 651;
soma_high_cut_off = 1012;
on_off_threshold = 831.5;

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
somaDict = Dict((df[_,1], vec(Array(df[_,2:4][xyz[1:2]]))) for _ in 1:size(df,1));
#display(somaDict);
#return
#somaDict = map!((x) -> (x[1], x[2][xyz[1:2]]), somaDict); # doesn't work

onSAC = [70244, 70243, 70242, 70241, 70240, 70239, 70238, 70237, 70236, 70235, 70234, 70233, 
		70232, 70231, 70230, 70229, 70228, 70227, 70225, 70224, 70223, 70222, 70221, 70220, 
		70219, 70218, 70217, 70216, 70215, 70214, 70213, 70212, 70211, 70209, 70208, 70207, 
		70206, 70205, 70204, 70203, 70202, 70201, 70200, 70199, 70198, 70197, 70196, 70195, 
		70194, 70193, 70192, 70191, 70189, 70188, 70187, 70186, 70185, 70184, 70183, 70182, 
		70181, 70180, 70179, 70178, 70176, 70174, 70172, 70171, 70170, 70169, 70168, 70164, 
		70163, 70162, 70161, 70029, 70028, 20204, 20196, 20062, 20044, 20032, 20030, 20025, 
		26007, 26009, 26010, 26011, 26012, 26013, 26014, 26015, 26016, 26017, 26139, 26143, 
		26169, 26174, 26180, 26183, 26184, 26186, 26197];

offSAC = [70167, 70158, 70156, 70155, 70154, 70152, 70151, 70150, 70149, 70148, 70147, 70146, 
		70145, 70141, 70138, 70137, 70134, 70133, 70131, 70130, 70129, 70128, 70127, 70126, 
		70125, 70124, 70123, 70122, 70121, 70120, 70119, 70118, 70117, 70116, 70115, 70114, 
		70113, 70112, 70111, 70110, 70109, 70108, 70106, 70102, 70100, 70099, 70096, 70095, 
		70093, 70090, 70089, 70088, 70087, 70086, 70085, 70084, 70083, 70082, 70081, 70080, 
		70079, 70077, 70076, 70068, 70066, 70050, 70048, 70035, 70034, 70033, 70032, 70031, 
		70030, 70027, 70026, 70025, 70024, 70023, 70016, 70014];

assert(length(onSAC) + length(offSAC) == length(somaDict));

# index 1 = ON, 2 = OFF 
onOffDict = Dict((_, 1) for _ in onSAC);
merge!(onOffDict, Dict((_, 2) for _ in offSAC));
assert(length(onOffDict) == length(somaDict));

h5fid = h5open(h5path, "r");
voldata = h5fid["/main"];
#% dims = size(data)  #(2048,5376,3456)
voldim = [size(voldata)...];
griddim = ceil(Int, (voldim[xyz[1:2]]-2) ./ blocksize);

#T = eltype(voldata);
grid_GCs = [Set{eltype(voldata)}() for onoff in 1:2, ii in 1:griddim[1], jj in 1:griddim[2]];

grid_denom = zeros(Int, 360, 2, griddim...);	# 360deg * on_off * ...
display(size(grid_denom))

@time begin
gc_num = Dict(0 => zeros(Int, 360, 2));

exceptions = [];

#for gridy = 1:100 #:griddim[2]
#	for gridx = 1:10 #:griddim[1]
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

		#for z = soma_low_cut_off+1:soma_high_cut_off-1	# has potential of cutting off GC voxels but we probably don't care
		for z = 2 : soma_high_cut_off-soma_low_cut_off	# has potential of cutting off GC voxels but we probably don't care
			onoff = z > on_off_threshold-soma_low_cut_off ? 2 : 1;
			gclist = grid_GCs[onoff, gridx, gridy];

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

							if onoff != onOffDict[cur]
								# should never happen
								coord = [[x, y] + grid_offset, z + soma_low_cut_off];
								warn("On off discrepency: $cur $onoff $(onOffDict[cur]) $coord")
								push!(exceptions, [cur, onoff, coord]);
							else

								#grid_denom(angle, gridx, gridy) += 1; # WTH: this creates a generic function....
								grid_denom[angle, onoff, gridx, gridy] += 1;

								# numerator
								contactinggcs = unique([
									 blockdata[x-1, y, z],  blockdata[x, y-1, z], blockdata[x, y, z-1],
									 blockdata[x+1, y, z],  blockdata[x, y+1, z], blockdata[x, y, z+1]]);
								for gc in contactinggcs
									get!(gc_num, gc, zeros(Int, 360, 2))[angle, onoff] += 1;
								end
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
allGCs = copy(grid_GCs[1]);
for gclist in grid_GCs
	union!(allGCs, gclist);
end
println("n GCs = ", length(allGCs));

end # @time
@time if true
gc_denom = Dict( (_, zeros(Int, 360, 2)) for _ in allGCs );
end # @time
@time if true
for _ in 1:length(grid_GCs[1,:])
	for onoff = 1:2
		for gc in grid_GCs[onoff,_]
			gc_denom[gc][:,onoff] += grid_denom[:,onoff,_];
		end
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
#gc_denom_mat = hcat([[elem...] for elem in gc_denom]...);
#gc_num_mat = hcat([[elem...] for elem in gc_num]...);
gc_denom_mat = (collect(keys(gc_denom)), collect(values(gc_denom)));
gc_num_mat = (collect(keys(gc_num)), collect(values(gc_num)));
#display(typeof(gc_denom_mat))

#file = matopen("gc_sac_contacts.mat", "w")
matwrite("gc_sac_contacts.mat", Dict(
	"gc_denom_keys" => gc_denom_mat[1],
	"gc_denom_vals" => gc_denom_mat[2],
	"gc_num_keys" => gc_num_mat[1],
	"gc_num_vals" => gc_num_mat[2],
	));
save("gc_sac_contacts.jld", Dict("gc_denom" => gc_denom, "gc_num" => gc_num, "exceptions" => exceptions));

return gc_num, gc_denom, exceptions

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
