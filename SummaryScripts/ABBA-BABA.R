## Function to compute ABBA-BABA statistic and z-score

ABBABABA = function(sample1) {
	
	n = sum(sample1[1,])
	p4 = sample1[,4]/n
	p7 = sample1[,7]/n
	p9 = sample1[,9]/n
	
	D = (p4-p7)/(p4+p7)
	D2 = (p9-p7)/(p9+p7)
	
	z = D/(2*sqrt(0.25/(n*p4+n*p7)))
	z2 = D2/(2*sqrt(0.25/(n*p9+n*p7)))

	return(cbind(D,z))
}