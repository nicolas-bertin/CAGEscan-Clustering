all: CAGEscan-Clustering.pod CAGEscan-Clustering.readme

CAGEscan-Clustering.pod: CAGEscan-Clustering.pl
	perldoc -u CAGEscan-Clustering.pl > CAGEscan-Clustering.pod

CAGEscan-Clustering.readme: CAGEscan-Clustering.pl
	perldoc -tT  CAGEscan-Clustering.pl > CAGEscan-Clustering.readme
