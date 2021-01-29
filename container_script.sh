docker pull nfcore/viralrecon:1.1.0
docker save nfcore/viralrecon:1.1.0 -o nfcore_viralrecon_1.1.0_docker.tar
singularity build nfcore_viralrecon_1.1.0.sif docker-archive://nfcore_viralrecon_1.1.0_docker.tar

docker build -t pangolin_nsc:latest - < pangolin_Dockerfile
docker save pangolin_nsc:latest -o pangolin_nsc_latest_docker.tar
singularity build pangolin_nsc_latest.sif docker-archive://pangolin_nsc_latest_docker.tar
