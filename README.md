# ap2

## getting started

1. Clone or download the repo.
2. Modify `config.yaml` to include your data.

## to get help

To print a description of every rule:

```bash
snakemake -l
```

To print a walkthrough of a google cloud run:

```bash
snakemake help
```

# local mode
```bash
snakemake --configfile config.yaml \
	  --use-singularity \
	  --directory /scratch/username/odir \
	  --keep-remote
```

# amarel mode
```bash
snakemake -j 50 --cluster-config amarel.json \
	  --cluster "sbatch --export=ALL --partition {cluster.partition} --nodes {cluster.n} --time {cluster.time} --ntasks {cluster.tasks} --mem {cluster.mem} --cpus-per-task={cluster.cpus}" \
	  --keep-remote \
	  --use-singularity \
	  --configfile config.yaml \
	  --directory /scratch/usename/odir
```

# GCP mode
```bash
# install prereqs
pip install kubernetes
gcloud components install kubectl

# set up cluster variables
CLUSTER_NAME=snk-cl2
NODES=8
ZONE=us-central1-a
REMOTE=GS
PREFIX=bucket_name
MACHINE_TYPE=n1-standard-2

# initialize cluster
gcloud container clusters create $CLUSTER_NAME \
	--num-nodes=$NODES \
	--scopes storage-rw \
	--machine-type=$MACHINE_TYPE \
	--zone $ZONE

# register cluster info
gcloud container clusters get-credentials $CLUSTER_NAME --zone $ZONE


snakemake --kubernetes --use-conda \
		--default-remote-provider $REMOTE \
		--default-remote-prefix $PREFIX \
		--latency-wait 300 \
		--jobs 8 \
		--verbose

# shut down your cluster
gcloud container clusters delete $CLUSTER_NAME --zone $ZONE
```
