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

# For any of the files in the dest_seqs folder (or its subfolders), we turn the dest
# seqs into CuspedOrbifold objects.
with open("./dest_seqs/Q_pi tilde covers, pi composite/CoversOfO20_5.txt", "r") as read_file:
    dest_seqs = read_file.read()


