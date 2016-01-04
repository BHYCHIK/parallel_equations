import random
N = 307199
for i in xrange(0, N):
	for j in xrange(0, N):
		if i == j:
			print random.uniform(-1000,1000), " ",
		else:
			print 0, " ",
	print random.uniform(-1000,1000)
		
