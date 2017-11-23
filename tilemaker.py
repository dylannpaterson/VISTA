import csv
import numpy as np

tiles = []
with open('wdb_query_20366_eso.csv') as csvfile:
	reader = csv.reader(csvfile)
	for row in reader:
		tiles.append(row)
print tiles[1]

tilear = np.asarray(tiles)
tilenames = list(set(tilear[:,1]))
tilenames.sort()

print tilenames

filenames = []



for tile in tilenames:
	tileid = np.where(tilear[:,1]==tile)[0]

	tilerange = tilear[tileid,2]

	print tile
	print tilear[tileid,1]

	seeingrange = tilerange.astype(float)
	
	seeingmin = np.min(seeingrange)

	filenames.append((tilear[tileid][np.where(tilear[tileid,2].astype(float)==seeingmin)[0],0][0]))


print len(set(tilear[:,0])), len(tilear[:,0])

outfile = open('test.txt', 'w')

for item in filenames:
  outfile.write("%s\n" % item)
  




