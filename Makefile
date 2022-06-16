##
# CSIDH
#
# @file
# @version 0.1

csidh_test:
	gcc -o sidh_test main.c fp.c mont_curve.c csidh.c -lgmp


# end
