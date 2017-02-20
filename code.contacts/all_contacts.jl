#Pkg.add("DataFrames")
#Pkg.add("MAT")
#Pkg.add("JLD")
using DataFrames
using MAT
using JLD

using HDF5;

function all_contacts(blocksize = [64; 64])

blocksize = vec(blocksize)

h5path = "/usr/people/smu/seungmount/research/Alex/retina/e2198_warped.h5";
#sacsomaCSV = "/data/smu/dev/e2198_Ca_imaging/data/sac_soma_centers_m2_warped.csv"
soma_low_cut_off = 600; #%90 #651;
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

bc = Dict(
	"bc1" => [60008, 60019, 60026, 60027, 60032, 60052, 60055, 60078, 
                   60079, 60099, 60105, 60109, 60110, 60111, 60114, 60118, 60129, 60132, 60139, 60142, 
                   60147, 60150, 60158, 60161, 60162, 60164, 60170, 60177, 60184, 60187, 60188, 60189, 
                   60194, 60195, 60196, 60203, 60204, 60212, 60213, 60216, 60218,],
	"bc2" => [60080, 60001, 60002, 60003, 60004, 60006, 60009, 60010, 
                   60011, 60012, 60013, 60021, 60022, 60023, 60025, 60029, 60037, 60038, 60039, 60040, 
                   60041, 60042, 60043, 60044, 60046, 60097, 60101, 60102, 60103, 60104, 60112, 60120, 
                   60124, 60130, 60133, 60135, 60138, 60140, 60141, 60149, 60157, 60159, 60167, 60169, 
                   60182, 60208, 60209, 60210, 60214, 60217, 60219, 60221, 60223, 60224, 60226,],
	"bc3a" => [60028, 60030, 60048, 60049, 60050, 60056, 60059, 60066, 
                   60068, 60072, 60075, 60076, 60085, 60088, 60092, 60093, 60108, 60115, 60123, 60134, 
                   60136, 60137, 60145, 60146, 60166, 60172, 60176, 60181,],
	"bc3b" => [60024, 60045, 60053, 60057, 60063, 60069, 60073, 60074, 
                   60077, 60082, 60084, 60087, 60090, 60094, 60096, 60098, 60106, 60107, 60122, 60131, 
                   60148, 60151, 60153, 60154, 60156, 60165, 60173, 60178, 60179, 60190, 60201, 60211, 
                   60215, 60228,],
	"bc4" => [60047, 60058, 60067, 60083, 60086, 60089, 60091, 60095, 
                   60100, 60116, 60117, 60119, 60121, 60125, 60127, 60143, 60144, 60160, 60163, 60168, 
                   60171, 60174, 60175, 60183, 60185, 60186, 60191, 60197, 60199, 60202, 60206, 60220, 
                   60222,],
	"bc5t" => [60543, 60452, 60366, 60371, 60387, 60395, 60399, 60403, 
                    60408, 60411, 60436, 60445, 60448, 60465, 60469, 60481, 60502, 60503, 60507, 60527, 
                    60533, 60535, 60538, 60054, 60155,],
	"bc5o" => [60364, 60390, 60405, 60406, 60420, 60432, 60434, 60440, 
                    60447, 60453, 60459, 60466, 60473, 60484, 60490, 60495, 60513, 60534, 60540, 60556, 
                    60559, 60608, 60612, 60018, 60061, 60064, 60071, 60532,],
	"bc5i" => [60360, 60033, 60374, 60380, 60383, 60386, 60388, 60389, 
                   60404, 60410, 60414, 60415, 60421, 60439, 60442, 60450, 60458, 60460, 60462, 60478, 
                   60488, 60491, 60497, 60498, 60504, 60505, 60510, 60514, 60519, 60522, 60523, 60528, 
                   60541, 60542, 60020, 60615, 60617, 60619, 60620, 60621, 60618,],
	"xbc" => [60355, 60449, 60379, 60413, 60430, 60455, 60493, 60501, 
                   60517, 60539, 60547, 60550,],
	"bc6" => [60489, 60422, 60356, 60060, 60423, 60394, 60425, 60398, 
                   60401, 60426, 60407, 60428, 60443, 60444, 60476, 60477, 60483, 60486, 60496, 60499, 
                   60512, 60516, 60526, 60548, 60549, 60552, 60554, 60561, 60611, 60614, 60467, 60358, 
                   60370, 60375, 60381, 60416, 60419, 60431, 60441, 60451, 60464, 60487, 60530, 60537, 
                   60036, 60363, 60409,],
	"bc7" => [60361, 60373, 60376, 60377, 60051, 60393, 60396, 60397, 
                   60412, 60418, 60429, 60435, 60437, 60446, 60454, 60470, 60480, 60485, 60492, 60508, 
                   60509, 60525, 60529, 60531, 60546, 60553, 60562, 60354, 60016,],
	"bc8/9" => [60433, 60482, 60500, 60578, 60368, 60402, 60417, 60438, 
                     60457, 60461,],
	"rbc" => [60463, 60536, 60357, 60359, 60365, 60369, 60372, 60378, 
                   60382, 60384, 60392, 60400, 60456, 60468, 60471, 60472, 60474, 60475, 60479, 60506, 
                   60515, 60518, 60520, 60521, 60544, 60551, 60555, 60557, 60558, 60560, 60563, 60564, 
                   60565, 60566, 60567, 60568, 60569, 60570, 60571, 60572, 60573, 60574, 60575, 60576, 
                   60577, 60580, 60581, 60582, 60583, 60584, 60585, 60586, 60587, 60588, 60589, 60590, 
                   60591, 60592, 60593, 60594, 60595, 60596, 60597, 60598, 60599, 60600, 60601, 60602, 
                   60603, 60604, 60605, 60606, 60607, 60609, 60610, 60613, 60017, 60031, 60427,],
);


bcDict = Dict((_, "bc1") for _ in bc["bc1"]);
for (name, cells) in bc
	merge!(bcDict, Dict((_, name) for _ in cells));
end

h5fid = h5open(h5path, "r");
voldata = h5fid["/main"];
#% dims = size(data)  #(2048,5376,3456)
voldim = [size(voldata)...];

griddim = ceil(Int, (voldim[xyz[1:2]]-2) ./ blocksize);

confusion_mat = Dict{UInt64, UInt}();

@time begin

exceptions = [];

#for gridy = 34:34 #:griddim[2]
#	for gridx = 34:34 #:griddim[1]
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
			
			for y = 2:jmax-1
				for x = 2:imax-1
					cur = blockdata[x, y, z];
					if cur == 0
						continue
					end
					thisKey = UInt64(cur) << 32;
                    # Caveat: This could over-count..especially for surface not aligned with axes, but actually
                    # corrects (but over corrects) the digitization error from using voxel count as a proxy for
                    # surface area
                    for neigbor in (blockdata[x+1, y, z],  blockdata[x, y+1, z],  blockdata[x, y, z+1])
                        if cur != neigbor && neigbor != 0
                            key = thisKey + neigbor;
							confusion_mat[key] = get(confusion_mat, key, 0) + 1;
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


#confusion_mat = (collect(keys(confusion_mat)), collect(values(confusion_mat)));
confusion_mat = hcat([[key >> 32, key << 32 >> 32, count] for (key, count) in confusion_mat]...);

display(typeof(confusion_mat))

#file = matopen("gc_bc_contacts.mat", "w")
matwrite("all_contacts.mat", Dict(
	"confusion_mat" => confusion_mat,
	));
save("all_contacts.jld", Dict("confusion_mat" => confusion_mat,
	"exceptions" => exceptions,));

return confusion_mat, exceptions

end # function

