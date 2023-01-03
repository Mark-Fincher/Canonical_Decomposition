"""
This is where we apply Sym to the orbifold census project.

1. Fix a Q pi tilde. Get the file with dests seqs of orbs in C_main which cover that Q pi tilde.
2. Turn those dest seqs into CuspedOrbifold objects.
3. Get all of their canonical triangulations.
4. Store each orbifold and isometry group in some way.
5. Determine how the isometry group permutes the lifts of the cusps of Q pi tilde.
6. Repeat this process for every other Q pi tilde.

With all the orbifold quotients, see which orbs we cannot tell apart just
with the singular locus.
"""

from Test import*

# For any of the files in the dest_seqs folder (or its subfolders), we turn the dest
# seqs into CuspedOrbifold objects.
with open("./dest_seqs/96seqs_sorted.txt", "r") as read_file:
    dest_seqs_string = read_file.read()

with open("48seqs.json", "r") as read_file:
    OrbDictionary = json.load(read_file)
    keyz = OrbDictionary.keys()
    OrbDictionary = {eval(k):OrbDictionary[k] for k in keyz}

ints_as_strings = ['0','1','2','3','4','5','6','7','8','9']

L = dest_seqs_string
dests = []
for i in range(len(L)):
	if L[i] == '[':
		dest = []
		getting_dest = True
	if L[i] == ']':
		dests.append(dest)
		getting_dest = False
		continue
	if L[i] in ints_as_strings:
		continue
	if L[i] == ',':
		continue
	if getting_dest == False:
		continue
	j = 1
	while L[i + j] in ints_as_strings:
		j += 1
	num = int(L[i+1:i+j])
	dest.append(num)

"""
failed = []
for dest in dests:
	print(' ')
	print(dest)
	orb = dest_to_orb(dest)
	if proto_canonize(orb):
		print('proto_canonize succeeded')
	else:
		failed.append(dests.index(dest))
if len(failed) > 0:
	print('failed for orbs with the following dest indices')
	print(failed)
else:
	print('proto canonize succeeded for all orbs')
"""

dests_failed = [107, 121, 180, 244, 245, 266, 337, 428, 429, 464, 485, 493, 523, 526, 530, 531, 538, 541, 556, 557, 564, 571, 572, 573, 587, 
588, 636, 642, 653, 723, 746, 777, 778, 783, 800, 802, 839, 840, 846, 850, 855, 857, 860, 861, 901, 906, 907, 918, 929, 930, 942, 943, 953, 
957, 958, 977, 984, 995, 997, 1000, 1001, 1002, 1038, 1040, 1266, 1267, 1268, 1269, 1270, 1274, 1277, 1296, 1308, 1311, 1312, 1313, 1315, 1320, 
1321, 1326, 1333, 1337, 1339, 1344, 1352, 1353, 1386, 1403, 1404, 1411, 1412, 1454, 1460, 1465, 1491, 1492, 1494, 1495, 1509, 1510, 1513, 1522, 
1523, 1524, 1573, 1586, 1593, 1624, 1625, 1693, 1770, 1772, 1774, 1829, 1830, 1834, 1838, 1839, 1840, 1841, 1842, 1866, 1933, 1934, 1942, 1950, 
1957, 1973, 1977, 1978, 1985, 1986, 1996, 2003, 2229, 2239, 2261, 2262, 2269, 2270, 2273, 2274, 2283, 2284, 2285, 2287, 2288, 2291, 2324, 2325, 
2326, 2327, 2332, 2333, 2341, 2343, 2351, 2352, 2366, 2368, 2370, 2371, 2372, 2373, 2374, 2376, 2378, 2381, 2382, 2383, 2403, 2432, 2433, 2463, 
2464, 2469, 2471, 2472, 2476, 2478, 2479, 2488, 2489, 2512, 2513, 2521, 2523, 2534, 2535, 2537, 2538, 2540, 2550, 2568, 2569, 2570, 2571, 2573, 
2574, 2587, 2589, 2597, 2599, 2631, 2632, 2633, 2635, 2640, 2641, 2642, 2643, 2644, 2650, 2653, 2672, 2674, 2697, 2698, 2717, 2720, 2721, 2735, 
2736, 2737, 2741, 2742, 2744, 2747, 2753, 2754, 2755, 2762, 2763, 2764, 2768, 2770, 2773, 2774, 2778, 2779, 2792, 2793, 2798, 2799, 2833, 2834, 
2835, 2836, 2837, 2838, 2839, 2840, 2841, 2842, 2843, 2844, 2845, 2846, 2847, 2848, 2852, 2853, 2860, 2863, 2864, 2865, 2866, 2867, 2868, 2870, 
2874, 2883, 2884, 2885, 2886, 2887, 2888, 2889, 2891, 2893, 2894, 2895, 2897, 2899, 2900, 2912, 2932, 2949, 2950, 2951, 2952, 2953, 2955, 2957, 
2958, 2959, 2961, 2963, 2976, 2997, 3017, 3019, 3022, 3031, 3036, 3047, 3049, 3057, 3059, 3068, 3069, 3076, 3077, 3094, 3095, 3109, 3110, 3111, 
3115, 3119, 3124, 3132, 3143, 3144, 3160, 3161, 3174, 3178, 3179, 3180, 3183, 3184, 3188, 3189, 3199, 3213, 3221, 3226, 3232, 3241, 3242, 3245, 
3252, 3268, 3269, 3282, 3283, 3288, 3325, 3330, 3331, 3336, 3349, 3392, 3400, 3403, 3404, 3408, 3409, 3410, 3428, 3433, 3434, 3436, 3451, 3476, 
3477, 3478, 3479, 3492, 3493, 3502, 3503, 3509, 3510, 3548, 3551, 3569, 3573, 3576, 3577, 3590, 3616, 3619, 3620, 3626, 3638, 3644, 3645, 3646, 
3647, 3648, 3650, 3653, 3682, 3684, 3698, 3703, 3757, 3761, 3801, 3803, 3809, 3810, 3819, 3820, 3834, 3835, 3837, 3842, 3846, 3847, 3876, 3877, 
3972, 4007, 4008, 4009, 4010, 4014, 4153]

"""
for failed_index in dests_failed:
	dest = dests[failed_index]
	orb = dest_to_orb(dest)
	if proto_canonize(orb):
		print("Proto canonize succeeded")
	else:
		print("Proto canonize failed")
	print(' ')
	print('---------------------------------------------------------------------')
	print(' ')
"""

dest = dests[107]
orb = dest_to_orb(dest)
print(proto_canonize(orb))
orb.info()
print(' ')
print_face_concavity(orb)
tet1 = orb.Tetrahedra[1]
tet0 = orb.Tetrahedra[0]
#print(orb.check_admissible(tet0))
#print(orb.special_four_to_four(tet1.Class[E02]))

"""
dest = [0, 0, 0, 1, 2, 3, 4, 0, 1, 5, 6, 2, 6, 7, 1, 4, 5, 1, 8, 3, 4, 9, 2, 10, 3, 2, 11, 12, 11, 13, 3, 7, 9, 4, 13, 13, 8, 14, 5, 15, 12, 12, 16, 5, 7, 6, 14, 17, 10, 18, 10, 6, 14, 8, 7, 8, 13, 11, 9, 19, 19, 20, 21, 9, 18, 10, 22, 20, 17, 23, 24, 11, 16, 25, 12, 24, 15, 26, 27, 14, 27, 28, 15, 16, 26, 15, 29, 26, 25, 16, 25, 30, 24, 31, 17, 27, 23, 17, 32, 18, 22, 22, 18, 33, 21, 34, 19, 21, 20, 19, 35, 23, 35, 36, 20, 37, 34, 21, 36, 38, 39, 40, 37, 22, 32, 41, 23, 42, 31, 24, 41, 43, 44, 43, 45, 25, 29, 46, 26, 47, 28, 27, 46, 48, 46, 29, 28, 49, 50, 30, 51, 28, 52, 53, 47, 29, 30, 50, 54, 44, 54, 55, 30, 45, 41, 32, 31, 56, 57, 48, 58, 31, 59, 60, 33, 32, 33, 61, 59, 39, 61, 33, 62, 40, 36, 35, 34, 63, 64, 38, 65, 34, 55, 54, 42, 35, 66, 67, 68, 36, 37, 58, 39, 55, 58, 37, 57, 67, 38, 64, 69, 64, 69, 70, 38, 68, 40, 39, 48, 71, 48, 57, 40, 50, 72, 62, 60, 41, 42, 51, 55, 59, 51, 42, 50, 62, 43, 44, 73, 57, 73, 56, 43, 60, 45, 74, 44, 75, 74, 45, 56, 58, 76, 77, 71, 46, 47, 78, 52, 52, 78, 47, 79, 77, 49, 68, 80, 76, 80, 81, 49, 51, 68, 49, 66, 53, 53, 52, 82, 69, 82, 79, 53, 83, 75, 63, 75, 54, 56, 73, 74, 72, 60, 59, 72, 81, 62, 72, 61, 80, 71, 71, 76, 61, 63, 75, 83, 66, 83, 83, 63, 65, 65, 84, 64, 82, 84, 65, 70, 79, 67, 66, 81, 74, 81, 80, 67, 73, 70, 69, 84, 78, 77, 76, 77, 70, 79, 82, 78, 85, 85, 85, 85, 84]
orb = dest_to_orb(dest)
print(proto_canonize(orb))
"""



"""
for key in OrbDictionary.keys():
	dest = OrbDictionary[key]
	orb = dest_to_orb(dest)
	if proto_canonize(orb):
		print('proto_canonize succeeded')
	else:
		raise Exception('proto_canonize failed on dest',dest)
"""	


# Covers of O(20,5).

"""

D = {}

D[(20, 5)] = [0, 1, 2, 3, 2, 4, 0, 2, 1, 0, 5, 1, 6, 3, 3, 0, 5, 7, 1, 8, 4, 2, 7, 9, 3, 6, 6, 6, 7, 5, 4, 10, 11, 9, 12, 4, 13, 12, 8, 5, 14, 15, 15, 7, 8, 16, 13, 13, 16, 8, 9, 15, 9, 11, 16, 11, 10, 17, 17, 14, 17, 10, 10, 12, 12, 13, 11, 18, 15, 14, 14, 19, 19, 18, 18, 16, 18, 19, 19, 17]
 
D[(40, 22)] = [1, 2, 2, 2, 0, 3, 3, 4, 3, 0, 0, 0, 2, 1, 1, 5, 5, 6, 7, 1, 4, 8, 9, 3, 9, 10, 4, 9, 8, 4, 11, 8, 7, 12, 5, 7, 6, 5, 13, 6, 13, 14, 6, 15, 12, 7, 14, 16, 11, 17, 8, 18, 10, 9, 17, 19, 17, 11, 10, 20, 21, 19, 22, 10, 23, 22, 18, 11, 14, 13, 12, 24, 25, 16, 26, 12, 27, 26, 15, 13, 28, 29, 29, 14, 15, 30, 27, 27, 30, 15, 16, 29, 16, 25, 30, 25, 31, 32, 32, 17, 18, 33, 23, 23, 33, 18, 19, 32, 19, 21, 33, 21, 20, 34, 34, 31, 34, 20, 20, 22, 22, 23, 21, 35, 24, 36, 36, 28, 36, 24, 24, 26, 26, 27, 25, 37, 29, 28, 28, 38, 38, 37, 37, 30, 32, 31, 31, 39, 39, 35, 35, 33, 35, 39, 39, 34, 37, 38, 38, 36]
 
D[(60, 130)] = [0, 1, 2, 3, 2, 4, 0, 2, 1, 0, 5, 1, 6, 3, 3, 0, 5, 7, 1, 8, 4, 2, 7, 9, 3, 6, 6, 6, 7, 5, 4, 10, 11, 9, 12, 4, 13, 14, 8, 5, 15, 16, 17, 7, 8, 18, 13, 13, 18, 8, 19, 16, 9, 11, 20, 11, 20, 21, 9, 17, 10, 22, 23, 15, 23, 24, 10, 12, 22, 10, 25, 14, 12, 26, 11, 27, 26, 12, 21, 28, 14, 13, 29, 30, 29, 19, 14, 31, 17, 32, 15, 33, 16, 15, 34, 35, 34, 36, 16, 37, 32, 17, 36, 38, 19, 29, 18, 29, 35, 30, 30, 18, 31, 39, 37, 19, 21, 20, 26, 26, 33, 27, 27, 20, 28, 38, 40, 21, 25, 41, 22, 34, 30, 35, 35, 22, 24, 23, 41, 32, 27, 33, 33, 23, 41, 25, 24, 42, 38, 28, 43, 24, 37, 44, 31, 25, 40, 45, 28, 40, 39, 31, 46, 39, 36, 34, 32, 47, 48, 49, 50, 36, 44, 37, 51, 49, 43, 52, 38, 50, 46, 51, 39, 53, 45, 40, 52, 54, 55, 47, 47, 41, 42, 56, 57, 55, 57, 50, 42, 43, 56, 42, 49, 44, 52, 43, 45, 52, 51, 46, 44, 51, 58, 54, 54, 45, 59, 53, 53, 46, 47, 55, 55, 48, 50, 57, 48, 57, 49, 48, 56, 56, 53, 59, 59, 59, 54, 58, 58, 58]
 
D[(60, 139)] = [0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 9, 0, 6, 10, 1, 11, 12, 13, 14, 1, 4, 2, 10, 15, 16, 17, 13, 2, 3, 18, 18, 8, 18, 3, 3, 13, 10, 6, 4, 19, 20, 21, 22, 4, 5, 23, 24, 16, 24, 7, 5, 9, 23, 5, 25, 21, 26, 22, 27, 6, 7, 24, 28, 12, 28, 25, 7, 27, 9, 8, 8, 29, 30, 31, 31, 10, 11, 32, 33, 26, 33, 34, 11, 14, 32, 11, 15, 31, 14, 35, 12, 36, 13, 12, 16, 37, 35, 14, 17, 38, 15, 39, 32, 20, 39, 15, 34, 17, 17, 16, 35, 40, 37, 29, 29, 18, 19, 41, 41, 30, 41, 19, 19, 22, 22, 26, 20, 42, 21, 20, 43, 44, 43, 27, 21, 45, 25, 28, 23, 46, 44, 40, 47, 23, 29, 37, 37, 24, 48, 45, 45, 25, 27, 43, 26, 49, 49, 47, 36, 28, 31, 30, 30, 50, 50, 42, 42, 32, 34, 33, 39, 51, 36, 52, 49, 33, 53, 38, 38, 34, 54, 55, 55, 35, 52, 36, 40, 55, 38, 53, 53, 54, 40, 44, 52, 39, 42, 50, 50, 41, 56, 57, 57, 43, 47, 49, 44, 57, 45, 48, 48, 56, 46, 58, 58, 48, 58, 46, 46, 47, 51, 59, 59, 53, 59, 51, 51, 52, 55, 54, 54, 59, 57, 56, 56, 58]
 
D[(60, 142)] = [0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 9, 0, 6, 10, 1, 11, 12, 13, 14, 1, 4, 2, 10, 15, 16, 17, 13, 2, 3, 18, 18, 8, 18, 3, 3, 13, 10, 6, 4, 19, 20, 21, 22, 4, 5, 23, 24, 16, 24, 7, 5, 9, 23, 5, 25, 21, 26, 27, 28, 6, 7, 24, 29, 12, 29, 25, 7, 28, 9, 8, 8, 30, 31, 32, 33, 10, 11, 34, 35, 26, 35, 22, 11, 14, 34, 11, 21, 32, 14, 36, 12, 37, 13, 12, 16, 38, 36, 14, 17, 39, 15, 40, 41, 20, 41, 37, 15, 33, 40, 15, 42, 17, 17, 16, 36, 43, 38, 30, 30, 18, 19, 44, 45, 31, 45, 39, 19, 22, 44, 19, 46, 27, 22, 35, 20, 47, 21, 20, 34, 48, 25, 29, 23, 46, 48, 43, 27, 23, 30, 38, 38, 24, 49, 50, 32, 25, 28, 51, 26, 52, 27, 26, 48, 53, 51, 28, 43, 50, 52, 42, 37, 29, 33, 54, 31, 55, 32, 31, 49, 56, 54, 33, 50, 36, 56, 53, 53, 34, 37, 41, 52, 35, 39, 45, 57, 54, 57, 46, 39, 42, 42, 52, 40, 58, 43, 48, 51, 40, 55, 47, 47, 41, 46, 57, 44, 49, 53, 56, 56, 44, 47, 55, 55, 45, 50, 49, 54, 59, 59, 58, 58, 51, 58, 59, 59, 57]
 
D[(60, 145)] = [0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 9, 0, 6, 10, 1, 11, 12, 13, 14, 1, 4, 2, 10, 15, 16, 17, 13, 2, 3, 18, 18, 8, 18, 3, 3, 13, 10, 6, 4, 19, 20, 21, 22, 4, 5, 23, 24, 16, 24, 7, 5, 9, 23, 5, 25, 21, 26, 27, 28, 6, 7, 24, 29, 12, 29, 25, 7, 28, 9, 8, 8, 30, 31, 32, 33, 10, 11, 34, 35, 26, 35, 36, 11, 14, 34, 11, 37, 32, 14, 38, 12, 39, 13, 12, 16, 40, 38, 14, 17, 41, 15, 42, 43, 20, 43, 28, 15, 33, 42, 15, 27, 17, 17, 16, 38, 37, 40, 30, 30, 18, 19, 44, 45, 31, 45, 46, 19, 22, 44, 19, 41, 27, 22, 47, 20, 48, 21, 20, 49, 50, 49, 39, 21, 51, 47, 22, 39, 29, 25, 29, 23, 46, 50, 37, 36, 23, 30, 40, 40, 24, 52, 33, 51, 25, 28, 43, 26, 47, 27, 26, 42, 53, 33, 52, 31, 54, 32, 31, 55, 56, 55, 51, 32, 38, 37, 50, 34, 42, 56, 53, 53, 34, 36, 35, 50, 57, 39, 49, 47, 35, 58, 41, 46, 36, 41, 58, 44, 55, 54, 48, 48, 43, 53, 56, 56, 44, 46, 45, 58, 52, 48, 54, 54, 45, 59, 57, 57, 49, 51, 55, 52, 59, 57, 59, 59, 58]
 
D[(60, 182)] = [0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 10, 0, 6, 11, 1, 12, 13, 14, 15, 1, 4, 2, 11, 16, 17, 18, 19, 2, 3, 20, 21, 8, 21, 10, 3, 19, 20, 3, 9, 14, 11, 6, 4, 22, 23, 24, 25, 4, 5, 26, 14, 17, 14, 13, 5, 10, 26, 5, 27, 24, 28, 29, 30, 6, 7, 19, 31, 13, 31, 32, 7, 30, 19, 7, 17, 9, 10, 21, 8, 20, 9, 8, 20, 21, 33, 34, 35, 11, 12, 36, 37, 28, 37, 25, 12, 15, 36, 12, 24, 34, 15, 27, 13, 38, 27, 15, 26, 39, 16, 40, 41, 23, 41, 30, 16, 35, 40, 16, 29, 18, 18, 17, 32, 42, 32, 31, 18, 43, 22, 44, 45, 33, 45, 39, 22, 25, 44, 22, 43, 29, 25, 37, 23, 46, 24, 23, 36, 47, 47, 42, 48, 26, 49, 50, 34, 27, 30, 41, 28, 51, 29, 28, 40, 52, 51, 48, 38, 31, 53, 35, 50, 32, 35, 53, 33, 54, 34, 33, 49, 55, 55, 52, 56, 36, 38, 57, 51, 37, 57, 38, 42, 50, 39, 45, 58, 49, 58, 43, 39, 48, 42, 47, 57, 40, 54, 56, 46, 41, 43, 58, 44, 53, 52, 55, 59, 44, 46, 59, 54, 45, 59, 46, 52, 57, 48, 51, 47, 56, 50, 49, 53, 59, 56, 54, 55, 58]
 
D[(80, 133)] = [0, 1, 2, 3, 2, 4, 0, 2, 1, 0, 5, 1, 6, 3, 3, 0, 5, 7, 1, 8, 4, 2, 7, 9, 3, 6, 6, 6, 7, 5, 4, 10, 11, 9, 12, 4, 13, 14, 8, 5, 15, 16, 17, 7, 8, 18, 13, 13, 18, 8, 19, 16, 9, 11, 20, 11, 20, 21, 9, 17, 10, 22, 23, 15, 23, 24, 10, 12, 22, 10, 25, 14, 12, 26, 11, 27, 26, 12, 21, 28, 14, 13, 29, 30, 29, 19, 14, 31, 17, 32, 15, 33, 16, 15, 34, 35, 34, 36, 16, 37, 32, 17, 36, 38, 19, 29, 18, 39, 35, 30, 40, 18, 41, 42, 37, 19, 21, 20, 26, 43, 33, 40, 27, 20, 44, 38, 45, 21, 25, 46, 22, 32, 30, 35, 47, 22, 24, 23, 46, 34, 27, 47, 33, 23, 46, 25, 24, 48, 37, 28, 41, 24, 38, 44, 31, 25, 49, 50, 51, 26, 47, 27, 30, 50, 28, 37, 52, 49, 52, 53, 28, 45, 54, 55, 50, 29, 31, 56, 38, 54, 56, 31, 57, 42, 36, 34, 32, 47, 40, 33, 35, 46, 50, 49, 54, 36, 39, 58, 48, 41, 48, 43, 39, 40, 58, 39, 59, 55, 42, 41, 53, 60, 53, 52, 42, 61, 43, 48, 62, 44, 62, 59, 43, 51, 45, 57, 44, 63, 57, 45, 56, 64, 51, 65, 49, 66, 65, 51, 55, 67, 63, 68, 66, 52, 69, 64, 70, 53, 55, 54, 65, 71, 60, 71, 68, 56, 72, 70, 61, 57, 59, 62, 58, 73, 71, 60, 74, 58, 75, 67, 67, 59, 68, 63, 60, 70, 61, 76, 72, 69, 76, 61, 64, 68, 66, 74, 63, 62, 64, 69, 76, 72, 77, 78, 78, 65, 74, 66, 71, 78, 67, 75, 75, 77, 70, 72, 69, 76, 73, 79, 79, 75, 79, 73, 73, 74, 78, 77, 77, 79]
 
D[(80, 164)] = [1, 0, 0, 2, 0, 1, 1, 3, 4, 5, 6, 0, 7, 8, 9, 1, 2, 10, 11, 7, 11, 6, 2, 6, 10, 2, 5, 5, 3, 12, 13, 4, 13, 14, 3, 9, 12, 3, 15, 8, 6, 11, 4, 16, 5, 4, 10, 17, 9, 18, 7, 19, 8, 7, 20, 21, 20, 22, 8, 23, 18, 9, 22, 24, 17, 21, 25, 10, 16, 26, 19, 11, 15, 27, 12, 28, 21, 17, 29, 12, 14, 13, 27, 30, 19, 31, 16, 13, 27, 15, 14, 32, 33, 24, 34, 14, 35, 34, 23, 15, 26, 16, 36, 26, 25, 37, 17, 25, 22, 20, 18, 38, 39, 40, 41, 18, 31, 19, 42, 40, 43, 44, 45, 20, 29, 46, 21, 45, 47, 48, 48, 22, 23, 49, 35, 43, 49, 23, 24, 48, 24, 33, 49, 39, 37, 25, 46, 50, 36, 42, 26, 51, 52, 53, 54, 27, 28, 55, 56, 35, 56, 41, 28, 29, 55, 28, 40, 53, 46, 29, 37, 57, 30, 58, 59, 33, 59, 45, 30, 54, 58, 30, 44, 31, 42, 36, 31, 60, 32, 61, 61, 52, 61, 32, 32, 34, 34, 35, 33, 62, 63, 64, 51, 36, 65, 50, 64, 37, 38, 66, 67, 47, 67, 57, 38, 41, 66, 38, 60, 44, 41, 56, 39, 68, 40, 39, 55, 58, 69, 70, 53, 42, 45, 59, 43, 56, 44, 43, 58, 71, 72, 54, 70, 46, 48, 47, 47, 73, 73, 71, 68, 49, 50, 65, 74, 65, 74, 51, 50, 70, 51, 74, 63, 63, 54, 72, 52, 75, 53, 52, 69, 76, 76, 62, 71, 55, 57, 67, 77, 72, 77, 60, 57, 64, 75, 68, 62, 59, 60, 77, 66, 69, 62, 76, 75, 61, 64, 63, 65, 78, 71, 73, 76, 66, 68, 75, 73, 67, 70, 69, 72, 79, 79, 78, 78, 74, 78, 79, 79, 77]
 
D[(80, 166)] = [1, 0, 0, 2, 0, 1, 1, 3, 4, 5, 6, 0, 7, 8, 9, 1, 2, 10, 11, 7, 11, 6, 2, 6, 10, 2, 5, 5, 3, 12, 13, 4, 13, 14, 3, 9, 12, 3, 15, 8, 6, 11, 4, 16, 5, 4, 10, 17, 9, 18, 7, 19, 8, 7, 20, 21, 20, 22, 8, 23, 18, 9, 22, 24, 17, 21, 25, 10, 16, 26, 19, 11, 15, 27, 12, 28, 21, 17, 29, 12, 14, 13, 27, 30, 19, 31, 16, 13, 27, 15, 14, 32, 33, 24, 34, 14, 35, 36, 23, 15, 26, 16, 37, 26, 25, 38, 17, 25, 22, 20, 18, 39, 40, 41, 36, 18, 31, 19, 42, 41, 43, 44, 45, 20, 29, 46, 21, 45, 47, 48, 49, 22, 23, 50, 35, 43, 50, 23, 41, 48, 24, 33, 51, 40, 51, 28, 24, 49, 38, 25, 46, 52, 37, 42, 26, 53, 54, 49, 55, 27, 28, 51, 56, 35, 56, 34, 28, 29, 46, 29, 38, 57, 30, 58, 59, 33, 59, 53, 30, 55, 58, 30, 60, 31, 42, 37, 31, 61, 32, 62, 63, 54, 63, 57, 32, 34, 62, 32, 39, 36, 34, 56, 33, 64, 36, 35, 40, 65, 66, 60, 53, 37, 67, 52, 44, 38, 39, 68, 62, 47, 68, 39, 57, 44, 41, 40, 50, 58, 69, 55, 48, 42, 45, 70, 43, 56, 44, 43, 67, 71, 70, 45, 52, 72, 73, 72, 72, 46, 49, 54, 47, 74, 48, 47, 69, 75, 75, 71, 76, 50, 74, 65, 65, 51, 52, 67, 70, 67, 53, 59, 66, 66, 55, 69, 54, 77, 57, 63, 68, 73, 60, 66, 58, 76, 77, 76, 64, 59, 78, 61, 61, 60, 61, 78, 78, 69, 65, 74, 74, 62, 64, 79, 77, 63, 79, 64, 71, 70, 71, 75, 79, 68, 72, 73, 73, 79, 76, 77, 75, 78]
 
D[(80, 168)] = [1, 0, 0, 2, 0, 1, 1, 3, 4, 5, 6, 0, 7, 8, 9, 1, 2, 10, 11, 7, 11, 6, 2, 6, 10, 2, 5, 5, 3, 12, 13, 4, 13, 14, 3, 9, 12, 3, 15, 8, 6, 11, 4, 16, 5, 4, 10, 17, 9, 18, 7, 19, 8, 7, 20, 21, 20, 22, 8, 23, 18, 9, 22, 24, 17, 21, 25, 10, 16, 26, 19, 11, 15, 27, 12, 28, 21, 17, 29, 12, 14, 13, 27, 30, 19, 31, 16, 13, 27, 15, 14, 32, 33, 24, 34, 14, 35, 36, 23, 15, 26, 16, 37, 26, 25, 38, 17, 25, 22, 20, 18, 39, 40, 41, 42, 18, 31, 19, 43, 41, 44, 34, 45, 20, 29, 46, 21, 45, 47, 48, 49, 22, 23, 50, 35, 44, 50, 23, 30, 48, 24, 33, 51, 40, 51, 45, 24, 49, 38, 25, 46, 52, 37, 43, 26, 53, 54, 55, 48, 27, 28, 56, 57, 35, 57, 58, 28, 29, 56, 28, 52, 55, 46, 29, 38, 59, 30, 60, 50, 33, 60, 30, 36, 31, 43, 37, 31, 61, 32, 62, 63, 54, 63, 39, 32, 34, 62, 32, 61, 36, 34, 44, 33, 64, 36, 35, 60, 65, 66, 42, 53, 37, 67, 52, 58, 38, 39, 63, 68, 47, 68, 61, 39, 42, 42, 66, 40, 69, 41, 40, 70, 60, 70, 53, 41, 71, 72, 71, 71, 43, 45, 51, 44, 57, 73, 49, 55, 46, 49, 73, 47, 74, 48, 47, 54, 75, 75, 64, 64, 50, 74, 76, 69, 51, 52, 67, 56, 67, 53, 70, 66, 66, 55, 54, 73, 77, 77, 65, 76, 56, 58, 57, 67, 76, 78, 59, 59, 58, 59, 78, 78, 73, 61, 68, 62, 72, 65, 77, 79, 62, 64, 75, 75, 63, 79, 69, 65, 70, 69, 79, 74, 68, 71, 72, 72, 79, 76, 74, 77, 78]
 
D[(80, 169)] = [1, 0, 0, 2, 0, 1, 1, 3, 4, 5, 6, 0, 7, 8, 9, 1, 2, 10, 11, 7, 11, 6, 2, 6, 10, 2, 5, 5, 3, 12, 13, 4, 13, 14, 3, 9, 12, 3, 15, 8, 6, 11, 4, 16, 5, 4, 10, 17, 9, 18, 7, 19, 8, 7, 20, 21, 20, 22, 8, 23, 18, 9, 22, 24, 17, 21, 25, 10, 16, 26, 19, 11, 15, 27, 12, 28, 21, 17, 29, 12, 14, 13, 27, 30, 19, 31, 16, 13, 27, 15, 14, 32, 33, 24, 34, 14, 35, 36, 23, 15, 26, 16, 37, 26, 25, 38, 17, 25, 22, 20, 18, 39, 40, 41, 42, 18, 31, 19, 43, 41, 44, 42, 45, 20, 29, 46, 21, 45, 47, 48, 49, 22, 23, 50, 35, 44, 50, 23, 51, 48, 24, 33, 52, 40, 52, 53, 24, 49, 38, 25, 46, 51, 37, 43, 26, 53, 54, 55, 55, 27, 28, 56, 57, 35, 57, 58, 28, 29, 56, 28, 30, 55, 46, 29, 38, 59, 30, 60, 56, 33, 60, 30, 58, 31, 43, 37, 31, 61, 32, 62, 63, 54, 63, 61, 32, 34, 62, 32, 59, 36, 34, 64, 33, 65, 64, 34, 53, 37, 36, 35, 66, 67, 66, 51, 36, 38, 39, 68, 68, 47, 68, 39, 39, 42, 42, 44, 40, 69, 41, 40, 70, 60, 70, 45, 41, 71, 72, 49, 71, 43, 45, 70, 44, 57, 73, 71, 48, 46, 49, 72, 47, 74, 48, 47, 73, 75, 51, 66, 50, 66, 75, 69, 67, 50, 53, 52, 64, 64, 74, 65, 69, 52, 55, 54, 54, 76, 76, 67, 65, 56, 58, 57, 60, 77, 78, 59, 61, 58, 59, 78, 62, 73, 61, 63, 78, 72, 67, 76, 75, 62, 65, 74, 76, 63, 69, 75, 74, 68, 79, 77, 77, 70, 71, 73, 72, 79, 77, 79, 79, 78]

# Just update these two lines, and D, for each case.
base_orb_str = '(20, 5)'
base_orb = dest_to_orb(D[(20,5)])

print("This file contains information about covers of", base_orb_str)
print(' ')
print('------------------------------------------------------------------')
print(' ')
for key in D.keys():
	print("Examining orbifold",key)
	orb = dest_to_orb(D[key])
	print(' ')
	print('Here is the starting triangulation')
	print(' ')
	orb.info()
	print(' ')
	print('The vertex classes of', key, 'are:')
	original_verts = [vertex for vertex in orb.Vertices]
	print(orb.Vertices)
	print('where')
	for vertex in orb.Vertices:
		print(vertex, '=', vertex.Corners)
	covering_maps = simplicial_maps_OP(orb,base_orb)
	print('There is/are', len(covering_maps), 'OP covering map(s) from', key, 'to', base_orb_str)
	map_0 = covering_maps[0]
	map_on_cusps = {}
	for tet in orb.Tetrahedra:
		for zero_subsimplex in ZeroSubsimplices:
			image_tet = map_0[tet][1]
			image_vertex = map_0[tet][0].image(zero_subsimplex)
			image_vertex_class = image_tet.Class[image_vertex]
			source_vertex_class = tet.Class[zero_subsimplex]
			map_on_cusps[source_vertex_class] = image_vertex_class
	print('The first covering map restricted to the vertex classes is:')
	print(map_on_cusps)
	print(' ')
	print("Attempting to make proto-canonical...")
	if proto_canonize(orb):
		print('Proto-canonize succeeded')
	else:
		print('Proto-canonize failed. Moving to next orb.')
		print(' ')
		print('---------------------------------------------------------------')
		print(' ')
		continue
	print(' ')
	print("Now doing canonize part 2...")
	canonical_retriangulation(orb)
	if orb.is_geometric is False:
		print("The proto-canonical triangulation was not the canonical decomp, so canonize part 2 found the canonical re-triangulation.")
	print(' ')
	print('The canonical triangulation is:')
	print(' ')
	orb.info()
	print(' ')
	print('Its vertex classes are')
	print(orb.Vertices)
	print('where')
	for vertex in orb.Vertices:
		if vertex in original_verts:
			print(vertex, 'is ideal and contains')
		else:
			print(vertex, 'is finite and contains')
		print(vertex.Corners)
	print(' ')
	print('We now print each OP self-isometry and how it acts on the vertex classes. The vertex classes were preserved')
	print('during the Pachner moves. e.g. if v0 was a cusp of the original triangulation, then v0 in this')
	print('triangulation is that same cusp.')
	for map_0 in simplicial_maps_OP(orb,orb):
		print(' ')
		print(map_0)
		map_on_classes = {}
		for tet in orb.Tetrahedra:
			for zero_subsimplex in ZeroSubsimplices:
				image_tet = map_0[tet][1]
				image_vertex = map_0[tet][0].image(zero_subsimplex)
				image_vertex_class = image_tet.Class[image_vertex]
				source_vertex_class = tet.Class[zero_subsimplex]
				map_on_classes[source_vertex_class] = image_vertex_class
		print('It acts on vertex classes as')
		print(map_on_classes)
	print(' ')
	print('-----------------------------------------------------------------------------------')
	print(' ')

"""




